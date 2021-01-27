#include <AMReX_Box.H>

#include <incflo.H>
#include <AMReX_NodalProjector.H>
#include <incflo_derive_K.H>

using namespace amrex;

void incflo::DiffFromExact (int /*lev*/, Geometry& lev_geom, Real time, Real dt,
                            MultiFab& error, int soln_comp, int err_comp) 
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
    // Decaying Taylor vortex
    } else if (2 == m_probtype) {

        for(MFIter mfi(error, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();

            // When we enter this routine, this holds the computed solution 
            Array4<Real> const& err = error.array(mfi);

            constexpr Real     pi =    3.1415926535897932;
            constexpr Real  twopi = 2.*3.1415926535897932;

            constexpr Real u0 = 1.;
            constexpr Real v0 = 1.;

            constexpr Real visc_coef = 0.001;

            Real omega = pi * pi * visc_coef;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {

                Real x = (i+0.5)*dx[0];
                Real y = (j+0.5)*dx[1];
                Real exact;
                if (err_comp == AMREX_SPACEDIM || err_comp == AMREX_SPACEDIM+1) {  // pressure 

                    Real t_p = time - 0.5*dt;
                    exact = -0.25 * ( std::cos(twopi*(x-u0*t_p)) + std::cos(twopi*(y-v0*t_p)) ) * std::exp(-4.*omega*t_p);
    
                } else if (err_comp == 0) { // u
                    exact =  u0 - std::cos(pi*(x-u0*time)) * std::sin(pi*(y-v0*time)) * std::exp(-2.*omega*time);
#if (AMREX_SPACEDIM == 3)
                    exact =  u0 - std::cos(pi*(x-u0*time)) * std::sin(pi*(y-v0*time)) * std::exp(-2.*omega*time);
#endif
    
                } else if (err_comp == 1) { // v
                    exact = v0 + std::sin(pi*(x-u0*time)) * std::cos(pi*(y-v0*time)) * std::exp(-2.*omega*time);
#if (AMREX_SPACEDIM == 3)
                    exact = v0 + std::sin(pi*(x-u0*time)) * std::cos(pi*(y-v0*time)) * std::exp(-2.*omega*time);
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
