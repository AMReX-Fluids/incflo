#include <incflo.H>

using namespace amrex;

//
// Analytic S(x,t) for each probtype: the default for compute_divu_source
//
void incflo::prob_divu_source (int lev, MultiFab& divu_source, Real time) const
{
    auto const& dx = geom[lev].CellSizeArray();
    auto const& problo = geom[lev].ProbLoArray();
    auto const& probhi = geom[lev].ProbHiArray();

    if (600 == m_probtype || 601 == m_probtype)
    {
        // S = f(t) prod_d cos(2 pi x_d / L_d), with f = 1 (600) or sin(2 pi t) (601).
        // In a periodic box, u = grad(phi) with lap(phi) = S at all times.
        constexpr Real twopi = Real(2.0) * Real(3.1415926535897932);
        Real const f = (600 == m_probtype) ? Real(1.0) : std::sin(twopi * time);
        for (MFIter mfi(divu_source, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& S = divu_source.array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                IntVect const iv(AMREX_D_DECL(i,j,k));
                Real s = f;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    Real const x = problo[d] + (Real(iv[d]) + Real(0.5)) * dx[d];
                    s *= std::cos(twopi * (x - problo[d]) / (probhi[d] - problo[d]));
                }
                S(i,j,k) = s;
            });
        }
    }
    else if (602 == m_probtype)
    {
        // Uniform S = 1
        divu_source.setVal(Real(1.0));
    }
    else
    {
        amrex::Abort("prob_divu_source: no S defined for this incflo.probtype");
    }
}
