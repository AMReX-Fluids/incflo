#include <incflo.H>

using namespace amrex;

Real
incflo::vol_wgt_sum (Vector<MultiFab*> const& mf_in, int icomp)
{
    Real  volwgtsum = 0.0;

    // Make a temporary copy so when we zero out levels we don't actually trash
    //  the original data
    Vector<MultiFab> mf(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(mf_in[lev]->boxArray(), mf_in[lev]->DistributionMap(), 1, 0);
        MultiFab::Copy(mf[lev], *mf_in[lev], icomp, 0, 1, 0);
    }

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const Real* dx = geom[lev].CellSize();

        // First, zero covered regions
        if (lev < finest_level)
        {
            BoxArray    baf;
            baf = mf_in[lev+1]->boxArray();
            baf.coarsen(ref_ratio[lev]);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                std::vector< std::pair<int,Box> > isects;
                for (MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                   auto const& fabarr = mf[lev].array(mfi);
                   int          ncomp = mf[lev].nComp();
                   baf.intersections(grids[lev][mfi.index()],isects);
                   for (const auto& is : isects)
                   {
                      amrex::ParallelFor(is.second, ncomp, [fabarr]
                      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                      {
                         fabarr(i,j,k,n) = 0.0;
                      });
                   } // isects
                } // mfi
            } // omp
        } // lev < finest

       // Use amrex::ReduceSum
       Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
#ifdef AMREX_USE_EB
       const EBFArrayBoxFactory* ebfact = &EBFactory(lev);
       auto const& vfrac = ebfact->getVolFrac();

       Real sm = amrex::ReduceSum(mf[lev], vfrac, 0, [vol]
       AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& mf_arr, Array4<Real const> const& vf_arr) -> Real
       {
           Real sum = 0.0;
           AMREX_LOOP_3D(bx, i, j, k,
           {
               sum += mf_arr(i,j,k) * vf_arr(i,j,k) * vol;
           });
           return sum;
       });
#else
       Real sm = amrex::ReduceSum(mf[lev], 0, [vol, dx]
       AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& mf_arr) -> Real
       {
           Real sum = 0.0;
           AMREX_LOOP_3D(bx, i, j, k,
           {
               sum += mf_arr(i,j,k) * vol;
           });
           return sum;
       });
#endif

        volwgtsum += sm;
    } // lev

    ParallelDescriptor::ReduceRealSum(volwgtsum);

    return volwgtsum;
}
