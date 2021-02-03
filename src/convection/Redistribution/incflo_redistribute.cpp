#ifdef AMREX_USE_EB
#include <Redistribution.H>

using namespace amrex;

void redistribution::redistribute_eb (Box const& bx, int ncomp,
                                      Array4<Real      > const& dUdt_out,
                                      Array4<Real const> const& dUdt_in,
                                      Array4<Real const> const& U_in,
                                      Array4<Real> const& scratch,
                                      AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                                   amrex::Array4<amrex::Real const> const& vmac,
                                                   amrex::Array4<amrex::Real const> const& wmac),
                                      Array4<EBCellFlag const> const& flag,
                                      AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                                   amrex::Array4<amrex::Real const> const& apy,
                                                   amrex::Array4<amrex::Real const> const& apz),
                                      amrex::Array4<amrex::Real const> const& vfrac,
                                      AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                                   amrex::Array4<amrex::Real const> const& fcy,
                                                   amrex::Array4<amrex::Real const> const& fcz),
                                      amrex::Array4<amrex::Real const> const& ccc,
                                      Geometry& lev_geom, Real dt, std::string redistribution_type)
{
    int redist_type;
    // redistribution_type = "NoRedist";          // no redistribution
    // redistribution_type = "FluxRedist";        // flux_redistribute
    // redistribution_type = "MergeRedistUpdate"; // merge_redistribute update
    // redistribution_type = "MergeRedistFull";   // merge_redistribute update
    // redistribution_type = "StateRedistUpdate"; // state_redistribute update
    // redistribution_type = "StateRedistFull";   // merge_redistribute full

#if (AMREX_SPACEDIM == 2)
    // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we  
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,1),4);
#else
    // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we  
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,1),8);
#endif
    itracker.setVal(0);

    amrex::Print() << "REDISTRIBUTION TYPE " << redistribution_type << std::endl;

    if (redistribution_type == "FluxRedist")
    {
        flux_redistribute_eb (bx, ncomp, dUdt_out, dUdt_in, scratch, flag, vfrac, lev_geom);

    } else if (redistribution_type == "MergeRedistUpdate") {
        Array4<int> itr = itracker.array();
        make_itracker(bx,
                      AMREX_D_DECL(apx, apy, apz), vfrac,
                      itr, lev_geom);

        merge_redistribute_update(bx, ncomp, dUdt_out, dUdt_in,
                                  AMREX_D_DECL(apx, apy, apz), vfrac,
                                  itr, lev_geom);

    } else if (redistribution_type == "MergeRedistFull") {
        Array4<int> itr = itracker.array();
        make_itracker(bx,
                      AMREX_D_DECL(apx, apy, apz), vfrac,
                      itr, lev_geom);

        merge_redistribute_full(bx, ncomp, dUdt_out, dUdt_in, U_in,
                                AMREX_D_DECL(apx, apy, apz), vfrac,
                                itr, lev_geom, dt);

    } else if (redistribution_type == "StateRedistUpdate") {
        state_redistribute_update(bx, ncomp, dUdt_out, dUdt_in, flag,
                                  AMREX_D_DECL(apx, apy, apz), vfrac,
                                  AMREX_D_DECL(fcx, fcy, fcz), ccc, lev_geom);

    } else if (redistribution_type == "StateRedistFull") {
        state_redistribute_full(bx, ncomp, dUdt_out, dUdt_in, U_in, flag,
                                AMREX_D_DECL(apx, apy, apz), vfrac,
                                AMREX_D_DECL(fcx, fcy, fcz), ccc, lev_geom, dt);

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
#endif
