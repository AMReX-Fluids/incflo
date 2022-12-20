/**
 * \file hydro_redistribution.cpp
 * \addtogroup Redistribution
 * @{
 *
 */

#include <hydro_redistribution.H>
#include <AMReX_EB_utils.H>

using namespace amrex;

void Redistribution::Apply ( Box const& bx, int ncomp,
                             Array4<Real      > const& dUdt_out,
                             Array4<Real      > const& dUdt_in,
                             Array4<Real const> const& U_in,
                             Array4<Real> const& scratch,
                             Array4<EBCellFlag const> const& flag_old,
			     Array4<EBCellFlag const> const& flag_new,
                             AMREX_D_DECL(Array4<Real const> const& apx_new,
                                          Array4<Real const> const& apy_new,
                                          Array4<Real const> const& apz_new),
			     Array4<amrex::Real const> const& vfrac_old,
                             Array4<amrex::Real const> const& vfrac_new,
                             AMREX_D_DECL(Array4<Real const> const& fcx,
                                          Array4<Real const> const& fcy,
                                          Array4<Real const> const& fcz),
                             Array4<Real const> const& ccc,
                             amrex::BCRec  const* d_bcrec_ptr,
                             Geometry const& lev_geom, Real dt,
                             std::string redistribution_type,
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
	dUdt_out(i,j,k,n) = 0.;
    });

    if (redistribution_type == "StateRedist")
    {
        Box const& bxg1 = grow(bx,1);
        Box const& bxg2 = grow(bx,2);
        Box const& bxg3 = grow(bx,3);
        Box const& bxg4 = grow(bx,4);

#if (AMREX_SPACEDIM == 2)
        // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we
        //    use the first component of this for the number of neighbors, so 4 comps
	// For Central Merging we include all surrounding cells, so 9 in 2D
        IArrayBox itracker(bxg4,9);
        // How many nbhds is a cell in
#else
        // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
        //    use the first component of this for the number of neighbors
        IArrayBox itracker(bxg4,8);
#endif
        FArrayBox nrs_fab(bxg3,1);
        FArrayBox alpha_fab(bxg3,2);

        // Total volume of all cells in my nbhd
        FArrayBox nbhd_vol_fab(bxg2,1);

        // Centroid of my nbhd
        FArrayBox cent_hat_fab     (bxg3,AMREX_SPACEDIM);

        Elixir eli_itr = itracker.elixir();
        Array4<int> itr = itracker.array();
        Array4<int const> itr_const = itracker.const_array();

        Elixir eli_nrs = nrs_fab.elixir();
        Array4<Real      > nrs       = nrs_fab.array();
        Array4<Real const> nrs_const = nrs_fab.const_array();

        Elixir eli_alpha = alpha_fab.elixir();
        Array4<Real      > alpha       = alpha_fab.array();
        Array4<Real const> alpha_const = alpha_fab.const_array();

        Elixir eli_nbf = nbhd_vol_fab.elixir();
        Array4<Real      > nbhd_vol       = nbhd_vol_fab.array();
        Array4<Real const> nbhd_vol_const = nbhd_vol_fab.const_array();

        Elixir eli_chf = cent_hat_fab.elixir();
        Array4<Real      > cent_hat       = cent_hat_fab.array();
        Array4<Real const> cent_hat_const = cent_hat_fab.const_array();

        Box domain_per_grown = lev_geom.Domain();
        AMREX_D_TERM(if (lev_geom.isPeriodic(0)) domain_per_grown.grow(0,1);,
                     if (lev_geom.isPeriodic(1)) domain_per_grown.grow(1,1);,
                     if (lev_geom.isPeriodic(2)) domain_per_grown.grow(2,1););

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

        amrex::ParallelFor(Box(scratch), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                const Real scale = (srd_update_scale) ? srd_update_scale(i,j,k) : Real(1.0);
                scratch(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n) / scale;

		if (i==9 && j == 8 ){
		    printf("scratch in: %e \n", scratch(i,j,k,n));
		}

            }
        );

        amrex::Print() << "Start itracker" << std::endl;
        
        MakeITracker(bx, flag_old, flag_new,
		     AMREX_D_DECL(apx_new, apy_new, apz_new), vfrac_new, 
                     itr, lev_geom, target_volfrac);
       
        amrex::Print() << "Start State Redistribution" << std::endl;

        MakeStateRedistUtils(bx, flag_old, flag_new, vfrac_new,
			     ccc, itr, nrs, alpha, nbhd_vol, cent_hat,
                             lev_geom, target_volfrac);
        
        StateRedistribute(bx, ncomp, dUdt_out, scratch, flag_old, flag_new, vfrac_old, vfrac_new,
                          AMREX_D_DECL(fcx, fcy, fcz), ccc,  d_bcrec_ptr,
                          itr_const, nrs_const, alpha_const, nbhd_vol_const,
                          cent_hat_const, lev_geom, srd_max_order);

        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                // Only update the values which actually changed -- this makes
                // the results insensitive to tiling -- otherwise cells that aren't
                // changed but are in a tile on which StateRedistribute gets called
                // will have precision-level changes due to adding/subtracting U_in
                // and multiplying/dividing by dt.   Here we test on whether (i,j,k)
                // has at least one neighbor and/or whether (i,j,k) is in the
                // neighborhood of another cell -- if either of those is true the
                // value may have changed
                
                if (i == 9 && j == 8){
                    amrex::Print() << "Pre dUdt_out" << IntVect(i,j) << dUdt_out(i,j,k,n) << std::endl;
                    amrex::Print() << "U_in: " << U_in(i,j,k,n) << std::endl;
                }

                if ( itr(i,j,k,0) > 0 || nrs(i,j,k) > 1.
		     || (!flag_old(i,j,k).isRegular() && flag_new(i,j,k).isRegular()) )
                {
                   const Real scale = (srd_update_scale) ? srd_update_scale(i,j,k) : Real(1.0);

                   if (!flag_old(i,j,k).isCovered()){
                       dUdt_out(i,j,k,n) = scale * (dUdt_out(i,j,k,n) - U_in(i,j,k,n)) / dt;
                   }
                }
                else
                {
                   dUdt_out(i,j,k,n) = dUdt_in(i,j,k,n);
                }

                if (i == 9 && j == 8)
                    amrex::Print() << "Post dUdt_out" << IntVect(i,j) << dUdt_out(i,j,k,n) << std::endl;
            }
        );

    } else if (redistribution_type == "NoRedist") {
        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dUdt_out(i,j,k,n) = dUdt_in(i,j,k,n);
            }
        );

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
                                     amrex::Array4<amrex::Real const> const& vfrac,
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
    std::string msg = "Redistribution::ApplyToInitialData: Shouldn't be here with redist type "+redistribution_type;
        amrex::Error(msg);
    }

    Box const& bxg2 = grow(bx,2);
    Box const& bxg3 = grow(bx,3);
    Box const& bxg4 = grow(bx,4);

#if (AMREX_SPACEDIM == 2)
    // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(bxg4,9);
#else
    // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(bxg4,8);
#endif
    FArrayBox nrs_fab(bxg3,1);
    FArrayBox alpha_fab(bxg3,2);

    // Total volume of all cells in my nbhd
    FArrayBox nbhd_vol_fab(bxg2,1);

    // Centroid of my nbhd
    FArrayBox cent_hat_fab  (bxg3,AMREX_SPACEDIM);

    Elixir eli_itr = itracker.elixir();
    Array4<int> itr = itracker.array();
    Array4<int const> itr_const = itracker.const_array();

    Elixir eli_nrs = nrs_fab.elixir();
    Array4<Real      > nrs       = nrs_fab.array();
    Array4<Real const> nrs_const = nrs_fab.const_array();

    Elixir eli_alpha = alpha_fab.elixir();
    Array4<Real      > alpha       = alpha_fab.array();
    Array4<Real const> alpha_const = alpha_fab.const_array();

    Elixir eli_nbf = nbhd_vol_fab.elixir();
    Array4<Real      > nbhd_vol       = nbhd_vol_fab.array();
    Array4<Real const> nbhd_vol_const = nbhd_vol_fab.const_array();

    Elixir eli_chf = cent_hat_fab.elixir();
    Array4<Real      > cent_hat       = cent_hat_fab.array();
    Array4<Real const> cent_hat_const = cent_hat_fab.const_array();

    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        U_out(i,j,k,n) = 0.;
    });

    MakeITracker(bx, flag, flag,
		 AMREX_D_DECL(apx, apy, apz), vfrac,
                 itr, lev_geom, target_volfrac);
    

    MakeStateRedistUtils(bx, flag, flag, vfrac, ccc, itr, nrs,
			 alpha, nbhd_vol, cent_hat, lev_geom, target_volfrac);
    

    StateRedistribute(bx, ncomp, U_out, U_in, flag, flag, vfrac, vfrac,
		      AMREX_D_DECL(fcx, fcy, fcz), ccc,  d_bcrec_ptr,
                      itr_const, nrs_const, alpha_const, nbhd_vol_const,
                      cent_hat_const, lev_geom, srd_max_order);
}
/** @} */
