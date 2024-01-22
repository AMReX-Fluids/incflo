#include <incflo.H>

using namespace amrex;

void incflo::prob_set_mixedBC_mask (Orientation ori, Box const& bx,
                                    Array4<Real> const& mask, int lev, Real /*time*/)
{
    if (1100 == m_probtype || 1101 == m_probtype || 1102 == m_probtype)
    {
        const int direction = ori.coordDir();
        int half_num_cells  = domain.length(direction) / 2;

        Orientation::Side side = ori.faceDir();
        if (side == Orientation::low) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (direction == 0) {
                    if (i <= half_num_cells) {
                        mask(i,j,k,0) = 0; // outflow on bottom
                    }
                }
                else if (direction == 1) {
                    if (j <= half_num_cells) {
                        mask(i,j,k,0) = 0;
                    }
                }
                else if (direction == 2) {
                    if (k <= half_num_cells) {
                        mask(i,j,k,0) = 0;
                    }
                }
            });
        } else {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (direction == 0) {
                    if (i > half_num_cells) {
                        mask(i,j,k,0) = 0; // outflow on top
                    }
                }
                else if (direction == 1) {
                    if (j > half_num_cells) {
                        mask(i,j,k,0) = 0;
                    }
                }
                else if (direction == 2) {
                    if (k > half_num_cells) {
                        mask(i,j,k,0) = 0;
                    }
                }
            });
        }
    }
    else
    {
        Abort("incflo::prob_set_mixedBC_mask: No masking function for probtype "
              +std::to_string(probtype));
    }
}
