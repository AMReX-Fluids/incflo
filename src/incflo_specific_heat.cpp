#include <incflo.H>

using namespace amrex;

void incflo::compute_cp (int /*lev*/, MFIter& /*mfi*/, FArrayBox& cp)
{
    // Get leveldata if desired, e.g.
    // Array4<Real const> const& rho   = density[lev]->array(mfi);

    Box const& bx = cp.box();
    Array4<Real> const& cp_a = cp.array();
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cp_a(i,j,k) = m_cp;
    });
}
