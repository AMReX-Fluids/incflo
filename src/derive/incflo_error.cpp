#include <AMReX_Box.H>

#include <incflo.H>
#include <AMReX_NodalProjector.H>
#include <derive_K.H>

using namespace amrex;

void incflo::DiffFromExact (int lev, Geometry& lev_geom, Real t, MultiFab& error, int soln_comp, int err_comp) 
{
    auto const& dx = lev_geom.CellSizeArray();

    // Taylor-Green vortices
    if (1 == m_probtype)
    {
        for(MFIter mfi(error, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();

            // When we enter this routine, this holds the computed solution 
            Array4<Real> const& err = error.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                constexpr Real  twopi = 2.*3.1415926535897932;
                constexpr Real fourpi = 4.*3.1415926535897932;

                Real x = (i+0.5)*dx[0];
                Real y = (j+0.5)*dx[1];
#if (AMREX_SPACEDIM == 3)
                Real z = (k+0.5)*dx[2];
#endif
                Real exact;
                if (err_comp == AMREX_SPACEDIM || err_comp == AMREX_SPACEDIM+1) {  // pressure 
                    exact = 0.25 * std::cos(fourpi*x) + 0.25 * std::cos(fourpi*y);
    
                } else if (err_comp == 0) { // u
                    exact =  std::sin(twopi*x) * std::cos(twopi*y);
#if (AMREX_SPACEDIM == 3)
                    exact *= std::cos(twopi*z);
#endif
    
                } else if (err_comp == 1) { // v
                    exact = -std::cos(twopi*x) * std::sin(twopi*y);
#if (AMREX_SPACEDIM == 3)
                    exact *= std::cos(twopi*z);
#endif

#if (AMREX_SPACEDIM == 3)
                } else if (err_comp == 2) { // w
                    exact = 0.;
#endif
                }

                err(i,j,k,soln_comp) -= exact;
            });
        }
    } else {
        amrex::Abort("Currently TGV is the only problem with an exact solution implemented");
    }
}
