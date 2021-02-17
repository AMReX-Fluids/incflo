#ifdef AMREX_USE_EB
#include <Redistribution.H>
#include <AMReX_EB_utils.H>

using namespace amrex;

void redistribution::redistribute_eb (Box const& bx, int ncomp,
                                      Array4<Real      > const& dUdt_out,
                                      Array4<Real      > const& dUdt_in,
                                      Array4<Real const> const& U_in,
                                      Array4<Real> const& scratch,
                                      Array4<EBCellFlag const> const& flag,
                                      AMREX_D_DECL(Array4<Real const> const& apx,
                                                   Array4<Real const> const& apy,
                                                   Array4<Real const> const& apz),
                                      Array4<amrex::Real const> const& vfrac,
                                      AMREX_D_DECL(Array4<Real const> const& fcx,
                                                   Array4<Real const> const& fcy,
                                                   Array4<Real const> const& fcz),
                                      Array4<Real const> const& ccc,
                                      Geometry& lev_geom, Real dt, std::string redistribution_type)
{
    // redistribution_type = "NoRedist";      // no redistribution
    // redistribution_type = "FluxRedist"     // flux_redistribute
    // redistribution_type = "MergeRedist";   // merge redistribute
    // redistribution_type = "StateRedist";   // state redistribute

#if (AMREX_SPACEDIM == 2)
    // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we  
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,1),4);
#else
    // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we  
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,1),8);
#endif

    Array4<int> itr = itracker.array();

    amrex::ParallelFor(Box(dUdt_out),ncomp, 
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            dUdt_out(i,j,k,n) = 0.;
        });

    if (redistribution_type == "FluxRedist")
    {
        amrex::ParallelFor(Box(scratch),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                scratch(i,j,k) = 1.;
            });

        int icomp = 0;
        apply_flux_redistribution (bx, dUdt_out, dUdt_in, scratch, icomp, ncomp, flag, vfrac, lev_geom);

    } else if (redistribution_type == "MergeRedist") {

        amrex::ParallelFor(Box(dUdt_in), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dUdt_in(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n);
            }
        );

        make_itracker(bx, AMREX_D_DECL(apx, apy, apz), vfrac, itr, lev_geom, "Merge");

        merge_redistribute(bx, ncomp, dUdt_out, dUdt_in, vfrac, itr, lev_geom);

        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dUdt_out(i,j,k,n) = (dUdt_out(i,j,k,n) - U_in(i,j,k,n)) / dt;
            }
        );

    } else if (redistribution_type == "StateRedist") {

        amrex::ParallelFor(Box(dUdt_in), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dUdt_in(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n);
            }
        );

        make_itracker(bx, AMREX_D_DECL(apx, apy, apz), vfrac, itr, lev_geom, "State");

        state_redistribute(bx, ncomp, dUdt_out, dUdt_in, flag, vfrac,
                           AMREX_D_DECL(fcx, fcy, fcz), ccc, itr, lev_geom);

        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dUdt_out(i,j,k,n) = (dUdt_out(i,j,k,n) - U_in(i,j,k,n)) / dt;
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

void redistribution::redistribute_initial_data (Box const& bx, int ncomp,
                                                Array4<Real      > const& U_inout,
                                                Array4<EBCellFlag const> const& flag,
                                                AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                                             amrex::Array4<amrex::Real const> const& apy,
                                                             amrex::Array4<amrex::Real const> const& apz),
                                                amrex::Array4<amrex::Real const> const& vfrac,
                                                AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                                             amrex::Array4<amrex::Real const> const& fcy,
                                                             amrex::Array4<amrex::Real const> const& fcz),
                                                amrex::Array4<amrex::Real const> const& ccc,
                                                Geometry& lev_geom, std::string redistribution_type)
{
    // redistribution_type = "MergeRedist";   // merge redistribute 
    // redistribution_type = "StateRedist";   // state redistribute

#if (AMREX_SPACEDIM == 2)
    // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we  
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,1),4);
#else
    // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we  
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,1),8);
#endif

    FArrayBox U_out(bx,ncomp);
    Array4<Real> uout_array = U_out.array(); 
    amrex::ParallelFor(Box(uout_array),ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        uout_array(i,j,k,n) = 0.;
    });

    Array4<int> itr = itracker.array();

    if (redistribution_type == "MergeRedist") {

        make_itracker(bx, AMREX_D_DECL(apx, apy, apz), vfrac, itr, lev_geom, "Merge");

        merge_redistribute(bx, ncomp, uout_array, U_inout, vfrac, itr, lev_geom);

    } else if (redistribution_type == "StateRedist") {

        make_itracker(bx, AMREX_D_DECL(apx, apy, apz), vfrac, itr, lev_geom, "State");

        state_redistribute(bx, ncomp, uout_array, U_inout, flag, vfrac,
                           AMREX_D_DECL(fcx, fcy, fcz), ccc, itr, lev_geom);

    } else {
       amrex::Error("Shouldnt be here with this redist type");
    }

    // Return the new (redistributed) data in the same U_in that came in
    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        U_inout(i,j,k,n) = uout_array(i,j,k,n);
    });
}
#endif
