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

// For moving SRD, fill newly uncovered cells in valid region of the MF with
// the value of it's merging neighbor
void
Redistribution::FillNewlyUncovered ( MultiFab& mf,
                                     EBFArrayBoxFactory const& ebfact_old,
                                     EBFArrayBoxFactory const& ebfact_new,
                                     MultiFab const& vel_eb,
                                     Geometry& geom,
                                     Real target_volfrac)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();

        EBCellFlagFab const& flagfab = ebfact_new.getMultiEBCellFlagFab()[mfi];
        if ( flagfab.getType(bx) == FabType::singlevalued)
        {
            AMREX_D_TERM(Array4<Real const> apx_old = ebfact_old.getAreaFrac()[0]->const_array(mfi);,
                         Array4<Real const> apy_old = ebfact_old.getAreaFrac()[1]->const_array(mfi);,
                         Array4<Real const> apz_old = ebfact_old.getAreaFrac()[2]->const_array(mfi););
            AMREX_D_TERM(Array4<Real const> apx_new = ebfact_new.getAreaFrac()[0]->const_array(mfi);,
                         Array4<Real const> apy_new = ebfact_new.getAreaFrac()[1]->const_array(mfi);,
                         Array4<Real const> apz_new = ebfact_new.getAreaFrac()[2]->const_array(mfi););
            Array4<Real const> vfrac_old = ebfact_old.getVolFrac().const_array(mfi);
            Array4<Real const> vfrac_new = ebfact_new.getVolFrac().const_array(mfi);
            Array4<Real const> vel_eb_arr= vel_eb.const_array(mfi);
            Array4<Real> U_in = mf.array(mfi);


// FIXME - how big does this box really need to be?
            // MakeITracker has 4 hard-coded into it, but here we would otherwise only need
            // 1 ghost cell
            Box const& gbx = grow(bx,4);

            IArrayBox itracker(gbx,itracker_comp,The_Async_Arena());
            Array4<int> itr = itracker.array();

            MakeITracker(bx, AMREX_D_DECL(apx_old, apy_old, apz_old), vfrac_old,
                             AMREX_D_DECL(apx_new, apy_new, apz_new), vfrac_new,
                         itr, geom, target_volfrac, vel_eb_arr);

            auto map = getCellMap();

            // Fill only valid region here. This will require FillPatch later...
            amrex::ParallelFor(Box(bx), mf.nComp(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                // Check to see if this cell was covered at time n, but uncovered at n+1
                if (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1. && vfrac_old(i,j,k) == 0.)
                {
                    for (int i_nbor = 1; i_nbor <= itr(i,j,k,0); i_nbor++)
                    {
                        int ioff = map[0][itr(i,j,k,i_nbor)];
                        int joff = map[1][itr(i,j,k,i_nbor)];
                        int koff = (AMREX_SPACEDIM < 3) ? 0 : map[2][itr(i,j,k,i_nbor)];

                        // Take the old value of my neighbor as my own
                        // NOTE this is only right for the case that the newly
                        // uncovered cell has only one other cell in it's neghborhood.
                        U_in(i,j,k,n) = U_in(i+ioff,j+joff,k+koff,n);

                        // FIXME -- correct fix of parallel OOB error here is that
                        // we check if we fall in the box...
                        amrex::Print() << "Cell  " << Dim3{i,j,k}
                                       << " newly uncovered, fill with value of neighbor at "
                                       << Dim3{i+ioff,j+joff,k+koff}
                                       <<": "<<U_in(i,j,k,n)<< std::endl;
                    }
                }
            });
        }
        else if ( !(flagfab.getType(bx) == FabType::regular || flagfab.getType(bx) == FabType::covered) )
        {
            Abort("Redistribution::FillNewlyUncovered(): Bad CellFlag type");
        }
    } //end mfiter
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

    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        out(i,j,k,n) = 0.;
    });

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

        amrex::Print() << "Start itracker" << std::endl;

        // FIXME - think about if this still needs v_eb and whether old or new...
        MakeITracker(bx, AMREX_D_DECL(apx_old, apy_old, apz_old), vfrac_old,
                         AMREX_D_DECL(apx_new, apy_new, apz_new), vfrac_new,
                     itr, lev_geom, target_volfrac, vel_eb_new);

        amrex::Print() << "Start State Redistribution" << std::endl;

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

            amrex::ParallelFor(Box(scratch), ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1. && vfrac_old(i,j,k) == 0.)
                {
                    // For SRD without slopes, shouldn't matter what's in here because
                    // it gets mult by V^n which is zero
                    // But, we need this to be consistent with how we define the update
                    // (out) below and in the application code
                    scratch(i,j,k,n) = U_in(i,j,k,n);
                    //scratch(i,j,k,n) = 0.0;
                }
                else if ((vfrac_old(i,j,k) > 0. && vfrac_old(i,j,k) < 1.0) ||
                         (vfrac_new(i,j,k) < 1. && vfrac_old(i,j,k) == 1.0) ||
                         (vfrac_old(i,j,k) == 1. && !flag_old(i,j,k).isRegular()))
                {
                    // Correct all cells that are cut at time n or become cut at time n+1

                    Real delta_vol = vfrac_new(i,j,k) - vfrac_old(i,j,k);
                    delta_vol /= ( dt * vfrac_old(i,j,k) );

                    // scratch(i,j,k,n) = 0.0;
                    // eb_add_divergence_from_flow(i,j,k,n,scratch,vel_eb,
                    //                             flag_old,vfrac_old,bnorm,barea,dxinv);

                    // Maybe we need to think of this correction as associated with the pred
                    // so for covering/uncovering we still need the full term and not half?
                    Real Ueb_dot_an =
                        AMREX_D_TERM(  vel_eb_old(i,j,k,0)*bnorm_old(i,j,k,0) * dxinv[0],
                                     + vel_eb_old(i,j,k,1)*bnorm_old(i,j,k,1) * dxinv[1],
                                     + vel_eb_old(i,j,k,2)*bnorm_old(i,j,k,2) * dxinv[2] );
                    Ueb_dot_an *= barea_old(i,j,k); ///vfrac_old(i,j,k);

                    if ( vel_eb_new && vfrac_new(i,j,k) > 0. ) {
                        Real Ueb_dot_an_new =
                            AMREX_D_TERM(  vel_eb_new(i,j,k,0)*bnorm_new(i,j,k,0) * dxinv[0],
                                         + vel_eb_new(i,j,k,1)*bnorm_new(i,j,k,1) * dxinv[1],
                                         + vel_eb_new(i,j,k,2)*bnorm_new(i,j,k,2) * dxinv[2] );
                        Ueb_dot_an_new *= barea_new(i,j,k); ///vfrac_new(i,j,k);

                        Ueb_dot_an = Real(0.5) * (Ueb_dot_an + Ueb_dot_an_new);

                    }

		    // This will undo volume scaling that happens later in forming q-hat
                    Ueb_dot_an /= vfrac_old(i,j,k);


                    Real delta_divU = delta_vol - Ueb_dot_an;
                    scratch(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n)
		    // FIXME - for now, comment delta-divU correction. This should be
		    // zero in 1D anyway.
                                                     // + U_in(i,j,k,n) * delta_vol
                                                     // - dt * Real(3.45) * Ueb_dot_an;
                                                     + dt * U_in(i,j,k,n) * delta_divU;

                    if ( (i==15 || i==16) && (j==11 || j==12) && n==0 ) {
                        Print()<<Dim3{i,j,k}<<" Uhat : "<<scratch(i,j,k,n) //<<std::endl;
                               // <<"Redist Pieces : "
			       <<" "<<U_in(i,j,k,n)
			       <<" "<<dUdt_in(i,j,k,n)*vfrac_old(i,j,k)<<" "
                               <<delta_divU<<" "
                               // <<Ueb_dot_an<<" "
                               // <<delta_vol
			       <<std::endl;
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
                            // amrex::Print() << "Cell  " << IntVect(i,j)
                            //                << " newly uncovered, correct neighbor at "
                            //                << Dim3{i+ioff,j+joff,k+koff} << std::endl;

                            Real delta_vol = vfrac_new(i,j,k) / vfrac_old(i+ioff,j+joff,k+koff);
                            // NOTE this correction is only right for the case that the newly
                            // uncovered cell has only one other cell in it's neghborhood.

                            // -eb_div * U_in *dt
                            // eb_add_divergence_from_flow(i,j,k,n,scratch,vel_eb,
                            //                             flag_old,vfrac_old,bnorm,barea,dxinv);

                            Real kappa_a = 0.;
                            if ( vel_eb_new) {
                                // Print()<<"Computing kappa "<<Dim3{i,j,k}<<n<<" from "
                                //        <<vel_eb_new(i,j,k,0)<<" "
                                //        <<bnorm_new(i,j,k,0)<<" "
                                //        <<barea_new(i,j,k)<<" "
                                //        <<out(i,j,k,n)<<" "
                                //        <<std::endl;
                                Real Ueb_dot_n =
                                    AMREX_D_TERM(  vel_eb_new(i,j,k,0)*bnorm_new(i,j,k,0) * dxinv[0],
                                                 + vel_eb_new(i,j,k,1)*bnorm_new(i,j,k,1) * dxinv[1],
                                                 + vel_eb_new(i,j,k,2)*bnorm_new(i,j,k,2) * dxinv[2] );

                                kappa_a = Real(0.5) * dt * Ueb_dot_n * barea_new(i,j,k) * U_in(i+ioff,j+joff,k+koff,n)
                                    / vfrac_old(i+ioff,j+joff,k+koff);

//                                 if ( j==8){
//                                     Print()<<"\nMSRD corrections...\n"
//                                            <<vel_eb(i,j,k,0)<<" "<<vel_eb(i,j,k,1)<<"\n"
//                                            <<bnorm(i,j,k,0)<<" "<<bnorm(i,j,k,1)<<"\n"
//                                            <<delta_vol<<" "<<kappa_a<<" "<<dUdt_in(i,j,k,n)
//                                            <<std::endl;

//FIXME new way to try to account for flux into newly uncovered cell.
				// This is needed for conservation
				//if (j==8) Print()<<"NU An+1 "<<dUdt_in(i,j,k,n)<<std::endl;
				scratch(i+ioff,j+joff,k+koff,n) += Real(0.5) * dt * dUdt_in(i,j,k,n) * vfrac_new(i,j,k)/vfrac_old(i+ioff,j+joff,k+koff);
			    }

                            // Print()<<"kappa "<<kappa_a<<" "<<delta_vol*U_in(i+ioff,j+joff,k+koff,n)<<std::endl;
                            scratch(i+ioff,j+joff,k+koff,n) += U_in(i+ioff,j+joff,k+koff,n) * delta_vol
				                               - kappa_a;
                        }
                    }
                }
            });
        }

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
                //fixme
                // if (i==8  && j==8){
                //     amrex::Print() << "Pre out" << Dim3{i,j,k} << out(i,j,k,n) << std::endl;
                //     amrex::Print() << "U_in: " << U_in(i,j,k,n) << std::endl;
                // }

                // FIXME - could probably make this logic more concise...
                if ( !( itr(i,j,k,0) > 0 || nrs(i,j,k) > 1.
                       || (vfrac_new(i,j,k) < 1. && vfrac_new(i,j,k) > 0.)
                       || (vfrac_old(i,j,k) < 1. && vfrac_new(i,j,k) == 1.) ) )
                {
                    // Only need to reset cells that didn't get SRD changes
                    out(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n);
                }

                //FIXME
                // if (i==0 && j==10)
                //     amrex::Print() << "Post out" << Dim3{i,j,k} << out(i,j,k,n) << std::endl;
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
