/**
 * \file hydro_redistribution.cpp
 * \addtogroup Redistribution
 * @{
 *
 */

#include <hydro_redistribution.H>
#include <AMReX_EB_utils.H>
#include <AMReX_EBMultiFabUtil_C.H>

using namespace amrex;

namespace {
    // For Normal Merging, we assume that in 2D a cell will need at most 3 neighbors to
    //   merge with. We use the first component of this for the number of neighbors, so
    //   4 comps needed.
    // For Central Merging, we include all surrounding cells, so in 2D, 9 comps needed.
    // For Moving EB, we have to allow for more then just Normal Merging (due to covering/
    //   uncovering reciprocity), so just allow for the max for now.
//
// FIXME --
// We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    constexpr int itracker_comp = (AMREX_SPACEDIM < 3 ) ? 9 : 8;
}

void Redistribution::Apply ( Box const& bx, int ncomp,
                             Array4<Real>       const& out,
                             Array4<Real>       const& dUdt_in,
                             Array4<Real const> const& U_in,
                             Array4<Real> const& scratch,
                             Array4<EBCellFlag const> const& flag,
                             AMREX_D_DECL(Array4<Real const> const& apx,
                                          Array4<Real const> const& apy,
                                          Array4<Real const> const& apz),
                             Array4<Real const> const& vfrac,
                             AMREX_D_DECL(Array4<Real const> const& fcx,
                                          Array4<Real const> const& fcy,
                                          Array4<Real const> const& fcz),
                             Array4<Real const> const& ccent,
                             BCRec  const* d_bcrec_ptr,
                             Geometry const& geom,
                             Real dt, std::string redistribution_type,
                             const int max_order,
                             Real target_volfrac,
                             Array4<Real const> const& update_scale)
{
    Apply(bx, ncomp, out, dUdt_in, U_in, scratch, flag, flag,
          AMREX_D_DECL(apx, apy, apz), vfrac,
          AMREX_D_DECL(apx, apy, apz), vfrac,
          AMREX_D_DECL(fcx, fcy, fcz), ccent,
          d_bcrec_ptr, geom, dt, redistribution_type,
          Array4<Real const> {}, // vel_eb_old
          Array4<Real const> {}, // bnorm_old
          Array4<Real const> {}, // barea_old, all not needed
          Array4<Real const> {}, // vel_eb
          Array4<Real const> {}, // bnorm
          Array4<Real const> {}, // barea, all not needed
          max_order, target_volfrac, update_scale);
}

void Redistribution::Apply ( Box const& bx, int ncomp,
                             Array4<Real      > const& out,
                             Array4<Real      > const& dUdt_in,
                             Array4<Real const> const& U_in,
                             Array4<Real> const& scratch,
                             Array4<EBCellFlag const> const& flag_old,
                             Array4<EBCellFlag const> const& flag_new,
                             AMREX_D_DECL(Array4<Real const> const& apx_old,
                                          Array4<Real const> const& apy_old,
                                          Array4<Real const> const& apz_old),
                             Array4<amrex::Real const> const& vfrac_old,
                             AMREX_D_DECL(Array4<Real const> const& apx_new,
                                          Array4<Real const> const& apy_new,
                                          Array4<Real const> const& apz_new),
                             Array4<amrex::Real const> const& vfrac_new,
                             AMREX_D_DECL(Array4<Real const> const& fcx,
                                          Array4<Real const> const& fcy,
                                          Array4<Real const> const& fcz),
                             Array4<Real const> const& ccc,
                             amrex::BCRec  const* d_bcrec_ptr,
                             Geometry const& lev_geom, Real dt,
                             std::string redistribution_type,
                             Array4<Real const> const& vel_eb_old,
                             Array4<Real const> const& bnorm_old,
                             Array4<Real const> const& barea_old,
                             Array4<Real const> const& vel_eb_new,
                             Array4<Real const> const& bnorm_new,
                             Array4<Real const> const& barea_new,
                             const int srd_max_order,
                             amrex::Real target_volfrac,
                             Array4<Real const> const& srd_update_scale)
{
    // redistribution_type = "NoRedist";       // no redistribution
    // redistribution_type = "FluxRedist"      // flux_redistribute
    // redistribution_type = "StateRedist";    // (weighted) state redistribute

//FIXME - For now, use the data in out
    // amrex::ParallelFor(bx,ncomp,
    // [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    // {
    //     out(i,j,k,n) = 0.;
    // });

    if (redistribution_type == "FluxRedist")
    {
        int icomp = 0;
        apply_flux_redistribution (bx, out, dUdt_in, scratch, icomp, ncomp,
                                   flag_old, vfrac_old, lev_geom);

    } else if (redistribution_type == "StateRedist") {

        Box const& bxg1 = grow(bx,1);
        Box const& bxg2 = grow(bx,2);
        Box const& bxg3 = grow(bx,3);
        Box const& bxg4 = grow(bx,4);

        // We use the first component of this for the number of neighbors, and later
        // components identify the neighbors (utilizing the CellMap)
        IArrayBox itracker(bxg4,itracker_comp,The_Async_Arena());

        FArrayBox nrs_fab(bxg3,1,The_Async_Arena());
        FArrayBox alpha_fab(bxg3,2,The_Async_Arena());

        // Total volume of all cells in my nbhd
        FArrayBox nbhd_vol_fab(bxg2,1,The_Async_Arena());

        // Centroid of my nbhd
        FArrayBox cent_hat_fab(bxg3,AMREX_SPACEDIM,The_Async_Arena());

        Array4<int> itr = itracker.array();
        Array4<int const> itr_const = itracker.const_array();

        Array4<Real      > nrs       = nrs_fab.array();
        Array4<Real const> nrs_const = nrs_fab.const_array();

        Array4<Real      > alpha       = alpha_fab.array();
        Array4<Real const> alpha_const = alpha_fab.const_array();

        Array4<Real      > nbhd_vol       = nbhd_vol_fab.array();
        Array4<Real const> nbhd_vol_const = nbhd_vol_fab.const_array();

        Array4<Real      > cent_hat       = cent_hat_fab.array();
        Array4<Real const> cent_hat_const = cent_hat_fab.const_array();

        Box domain_per_grown = lev_geom.Domain();
        AMREX_D_TERM(if (lev_geom.isPeriodic(0)) domain_per_grown.grow(0,1);,
                     if (lev_geom.isPeriodic(1)) domain_per_grown.grow(1,1);,
                     if (lev_geom.isPeriodic(2)) domain_per_grown.grow(2,1););

        if ( dUdt_in )
        {
            // At any external Dirichlet domain boundaries we need to set dUdt_in to 0
            //    in the cells just outside the domain because those values will be used
            //    in the slope computation in state redistribution.  We assume here that
            //    the ext_dir values of U_in itself have already been set.
            if (!domain_per_grown.contains(bxg1))
                amrex::ParallelFor(bxg1,ncomp,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (!domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k))))
                        dUdt_in(i,j,k,n) = 0.;
                });
        }


        // FIXME - think about if this still needs v_eb and whether old or new...
        MakeITracker(bx, AMREX_D_DECL(apx_old, apy_old, apz_old), vfrac_old,
                         AMREX_D_DECL(apx_new, apy_new, apz_new), vfrac_new,
                     itr, lev_geom, target_volfrac, vel_eb_new);


        MakeStateRedistUtils(bx, vfrac_old, vfrac_new, ccc, itr, nrs, alpha, nbhd_vol, cent_hat,
                             lev_geom, target_volfrac);

        if ( !vel_eb_old )
        {
            //
            // SRD with stationary EB
            //
            if ( dUdt_in )
            {
                // We're working with an update
                amrex::ParallelFor(Box(scratch), ncomp,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    const Real scale = (srd_update_scale) ? srd_update_scale(i,j,k) : Real(1.0);
                    scratch(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n) / scale;
                });
            }
            else
            {
                // We're doing a whole state
                amrex::ParallelFor(Box(scratch), ncomp,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    // FIXME? for this case I think we could do away with scratch and
                    // just use U_in
                    // Also, what about scale?
                    scratch(i,j,k,n) = U_in(i,j,k,n);
                });
            }
        }
        else
        {
            //
            // Moving SRD corrections
            //
            // Here, delta-divU is the difference between reasonable divU values
            // that we could pass into the MAC
            //
            const GpuArray<Real,AMREX_SPACEDIM> dxinv = lev_geom.InvCellSizeArray();
            auto map = getCellMap();

            // FIXME - for now, don't allow scaling with MSRD.
            AMREX_ALWAYS_ASSERT(!srd_update_scale);

	    Real eps = 1.e-14;
	    
            amrex::ParallelFor(Box(scratch), ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1. && vfrac_old(i,j,k) == 0.)
                {
                    // Newly Uncovered cells:
                    // For SRD without slopes, it shouldn't matter what's in here because
                    // it gets mult by V^n which is zero
                    scratch(i,j,k,n) = U_in(i,j,k,n);
                }
                else if ( (vfrac_old(i,j,k) > 0. && vfrac_old(i,j,k) < 1.0) ||
                         (vfrac_old(i,j,k) == 1. &&
                          (!flag_old(i,j,k).isRegular() || !flag_new(i,j,k).isRegular()) ))
                {
                    // Correct all cells that are cut at time n or become cut at time n+1
                    Real delta_divU = 0.0;
                    Real Ueb_dot_an = 0.0;
                    Real delta_vol = (vfrac_new(i,j,k) - vfrac_old(i,j,k))/dt;

                    if (!flag_old(i,j,k).isRegular() && !flag_new(i,j,k).isCovered()
                        && delta_vol > eps ) // we need this correction for NU but don't for NC...
                    {
                        Ueb_dot_an =
                            AMREX_D_TERM(  vel_eb_old(i,j,k,0)*bnorm_old(i,j,k,0) * dxinv[0],
                                           + vel_eb_old(i,j,k,1)*bnorm_old(i,j,k,1) * dxinv[1],
                                           + vel_eb_old(i,j,k,2)*bnorm_old(i,j,k,2) * dxinv[2] );
                        Ueb_dot_an *= barea_old(i,j,k);

                        delta_divU = (delta_vol - Ueb_dot_an) * U_in(i,j,k,n);
                    }

                    // For the Corrector step
                    if (!flag_new(i,j,k).isCovered() && delta_vol > eps && vel_eb_new)
                    {
                        Real Ueb_dot_an_new =
                            AMREX_D_TERM(  vel_eb_new(i,j,k,0)*bnorm_new(i,j,k,0) * dxinv[0],
                                         + vel_eb_new(i,j,k,1)*bnorm_new(i,j,k,1) * dxinv[1],
                                         + vel_eb_new(i,j,k,2)*bnorm_new(i,j,k,2) * dxinv[2] );
                        Ueb_dot_an_new *= barea_new(i,j,k);

                        if ( flag_old(i,j,k).isRegular() ){
                            delta_divU = 0.5 * (U_in(i,j,k,n) + out(i,j,k,n)) * delta_vol
                                - 0.5 * out(i,j,k,n) * Ueb_dot_an_new;

                             Abort("This not yet tested...");
                        } else {
                            delta_divU = Real(0.5) * (delta_divU
                                                      + out(i,j,k,n) * (delta_vol - Ueb_dot_an_new));
                        }
                    }

                    // This will undo volume scaling that happens later in forming q-hat
                    delta_divU /= vfrac_old(i,j,k);

                    scratch(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n)
                        + dt * delta_divU;

                    if (j==8 && (i==9 || i==8) ) {
                        //if (j==12 && (i==15 || i==16) ) {
                        Print()<<"NOT REGULAR! "<<i<<std::endl;
			Print()<<"DELTA_DIVU "<<i<<": "<<delta_vol<<" "<<Ueb_dot_an
                               <<" "<<delta_divU<<std::endl;
                     Print()<<U_in(i,j,k,n)<<std::endl;
                    }

                }
                else
                {
                    scratch(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n);
                }
            });

            // FIXME need to think about how big we really need this box to be
            amrex::ParallelFor(Box(scratch), ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                // Check to see if this cell was covered at time n.
                // If covered, add my vol correction to the cells in my nbhd
                // Need this in a separate loop because otherwise we can't be guaranteed
                // that the nb has been initialized.
                if ( vfrac_old(i,j,k) == 0.0 ) {
                    for (int i_nbor = 1; i_nbor <= itr(i,j,k,0); i_nbor++)
                    {
                        int ioff = map[0][itr(i,j,k,i_nbor)];
                        int joff = map[1][itr(i,j,k,i_nbor)];
                        int koff = (AMREX_SPACEDIM < 3) ? 0 : map[2][itr(i,j,k,i_nbor)];

                        if ( Box(scratch).contains(Dim3{i+ioff,j+joff,k+koff}) )
                        {
                            Real delta_vol = vfrac_new(i,j,k) / dt;
                            Real delta_divU = delta_vol * U_in(i+ioff,j+joff,k+koff,n);
			    Real dV = vfrac_new(i+ioff,j+joff,k+koff)-vfrac_old(i+ioff,j+joff,k+koff);
			    
                            // NOTE this correction is only right for the case that the newly
                            // uncovered cell has only one other cell in it's neghborhood.
			    if (!flag_old(i+ioff,j+joff,k+koff).isRegular() && !flag_new(i+ioff,j+joff,k+koff).isCovered()
				&& dV < eps ) // we need this correction for NU but don't for NC...
			    {
				Real Ueb_dot_an =
				    AMREX_D_TERM(  vel_eb_old(i+ioff,j+joff,k+koff,0)*bnorm_old(i+ioff,j+joff,k+koff,0) * dxinv[0],
						   + vel_eb_old(i+ioff,j+joff,k+koff,1)*bnorm_old(i+ioff,j+joff,k+koff,1) * dxinv[1],
						   + vel_eb_old(i+ioff,j+joff,k+koff,2)*bnorm_old(i+ioff,j+joff,k+koff,2) * dxinv[2] );
				Ueb_dot_an *= barea_old(i+ioff,j+joff,k+koff);

				delta_divU = (delta_vol - Ueb_dot_an) * U_in(i+ioff,j+joff,k+koff,n);
			    }


			    if ( vel_eb_new) {
                                Real Ueb_dot_n =
                                    AMREX_D_TERM(  vel_eb_new(i,j,k,0)*bnorm_new(i,j,k,0) * dxinv[0],
                                                 + vel_eb_new(i,j,k,1)*bnorm_new(i,j,k,1) * dxinv[1],
                                                 + vel_eb_new(i,j,k,2)*bnorm_new(i,j,k,2) * dxinv[2] );

                                Ueb_dot_n *= barea_new(i,j,k);

				delta_divU = 0.5*(delta_divU + 
						  out(i,j,k,n) * (delta_vol - Ueb_dot_n));
				
                                // Account for flux into newly uncovered cell (needed for conservation)
                                scratch(i+ioff,j+joff,k+koff,n) += Real(0.5) * dt * dUdt_in(i,j,k,n)
				    * vfrac_new(i,j,k)/vfrac_old(i+ioff,j+joff,k+koff);
                            }

                            scratch(i+ioff,j+joff,k+koff,n) += dt*delta_divU/vfrac_old(i+ioff,j+joff,k+koff);
                        }
                    }
                }
            });
        }

        //FIXME - For now, use the data in out, need to zero here
        amrex::ParallelFor(bx,ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            out(i,j,k,n) = 0.;
        });

        StateRedistribute(bx, ncomp, out, scratch, flag_new, vfrac_old, vfrac_new,
                          AMREX_D_DECL(fcx, fcy, fcz), ccc,  d_bcrec_ptr,
                          itr_const, nrs_const, alpha_const, nbhd_vol_const,
                          cent_hat_const, lev_geom, srd_max_order);

        //
        // Only update the values which actually changed -- this makes
        // the results insensitive to tiling -- otherwise cells that aren't
        // changed but are in a tile on which StateRedistribute gets called
        // will have precision-level changes due to adding/subtracting U_in
        // and multiplying/dividing by dt.   Here we test on whether (i,j,k)
        // has at least one neighbor and/or whether (i,j,k) is in the
        // neighborhood of another cell -- if either of those is true the
        // value may have changed
        //
        if ( !vel_eb_old )
        {
            //
            // SRD with stationary EB
            //
            if ( dUdt_in )
            {
                // Pass out an update
                amrex::ParallelFor(bx, ncomp,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (itr(i,j,k,0) > 0 || nrs(i,j,k) > 1.)
                    {
                        const Real scale = (srd_update_scale) ? srd_update_scale(i,j,k) : Real(1.0);

                        out(i,j,k,n) = scale * (out(i,j,k,n) - U_in(i,j,k,n)) / dt;

                    }
                    else
                    {
                        out(i,j,k,n) = dUdt_in(i,j,k,n);
                    }
                });
            }
            else
            {
                // Want to pass out the whole state, so we only need to reset cells that
                // didn't get SRD changes.
                amrex::ParallelFor(bx, ncomp,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if ( !(itr(i,j,k,0) > 0 || nrs(i,j,k) > 1.) )
                    {
                        out(i,j,k,n) = U_in(i,j,k,n);
                    }
                });
            }
        }
        else
        {
            //
            // MSRD - pass out the full redistributed state.
            //
            amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                // FIXME - could probably make this logic more concise...
                if ( !( itr(i,j,k,0) > 0 || nrs(i,j,k) > 1.
                       || (vfrac_new(i,j,k) < 1. && vfrac_new(i,j,k) > 0.)
                       || (vfrac_old(i,j,k) < 1. && vfrac_new(i,j,k) == 1.) ) )
                {
                    // Only need to reset cells that didn't get SRD changes
                    out(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n);
                }
        });
        }
    } else if (redistribution_type == "NoRedist") {
        Print()<<"No redistribution..."<<std::endl;

        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            out(i,j,k,n) = dUdt_in(i,j,k,n);
        });

    } else {
       amrex::Error("Not a legit redist_type");
    }
}

void
Redistribution::ApplyToInitialData ( Box const& bx, int ncomp,
                                     Array4<Real      > const& U_out,
                                     Array4<Real      > const& U_in,
                                     Array4<EBCellFlag const> const& flag,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                                  amrex::Array4<amrex::Real const> const& apy,
                                                  amrex::Array4<amrex::Real const> const& apz),
// FIXME - do we really need to both vfrac here?
                                     amrex::Array4<amrex::Real const> const& vfrac_old,
                                     amrex::Array4<amrex::Real const> const& vfrac_new,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                                  amrex::Array4<amrex::Real const> const& fcy,
                                                  amrex::Array4<amrex::Real const> const& fcz),
                                     amrex::Array4<amrex::Real const> const& ccc,
                                     amrex::BCRec  const* d_bcrec_ptr,
                                     Geometry& lev_geom, std::string redistribution_type,
                                     const int srd_max_order,
                                     amrex::Real target_volfrac)
{
    if (redistribution_type != "StateRedist") {
        std::string msg = "Redistribution::ApplyToInitialData: Shouldn't be here with redist type "
            +redistribution_type;
        amrex::Error(msg);
    }

    Box const& bxg2 = grow(bx,2);
    Box const& bxg3 = grow(bx,3);
    Box const& bxg4 = grow(bx,4);

// FIXME itracker comp should allow letting go of this #if
#if (AMREX_SPACEDIM == 2)
    // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(bxg4,itracker_comp,The_Async_Arena());
#else
    // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(bxg4,itracker_comp,The_Async_Arena());
#endif
    FArrayBox nrs_fab(bxg3,1,The_Async_Arena());
    FArrayBox alpha_fab(bxg3,2,The_Async_Arena());

    // Total volume of all cells in my nbhd
    FArrayBox nbhd_vol_fab(bxg2,1,The_Async_Arena());

    // Centroid of my nbhd
    FArrayBox cent_hat_fab(bxg3,AMREX_SPACEDIM,The_Async_Arena());

    Array4<int> itr = itracker.array();
    Array4<int const> itr_const = itracker.const_array();

    Array4<Real      > nrs       = nrs_fab.array();
    Array4<Real const> nrs_const = nrs_fab.const_array();

    Array4<Real      > alpha       = alpha_fab.array();
    Array4<Real const> alpha_const = alpha_fab.const_array();

    Array4<Real      > nbhd_vol       = nbhd_vol_fab.array();
    Array4<Real const> nbhd_vol_const = nbhd_vol_fab.const_array();

    Array4<Real      > cent_hat       = cent_hat_fab.array();
    Array4<Real const> cent_hat_const = cent_hat_fab.const_array();

    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        U_out(i,j,k,n) = 0.;
    });

    MakeITracker(bx, AMREX_D_DECL(apx, apy, apz), vfrac_old,
                     AMREX_D_DECL(apx, apy, apz), vfrac_new,
                 itr, lev_geom, target_volfrac);


    MakeStateRedistUtils(bx, vfrac_old, vfrac_new, ccc, itr, nrs, alpha, nbhd_vol, cent_hat,
                         lev_geom, target_volfrac);


    StateRedistribute(bx, ncomp, U_out, U_in, flag, vfrac_old, vfrac_new,
                      AMREX_D_DECL(fcx, fcy, fcz), ccc,  d_bcrec_ptr,
                      itr_const, nrs_const, alpha_const, nbhd_vol_const,
                      cent_hat_const, lev_geom, srd_max_order);
}
/** @} */
