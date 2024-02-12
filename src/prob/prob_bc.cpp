#include <incflo.H>

using namespace amrex;

// This function is used to prepare both the nodal projection and advection for
// mixed BCs (i.e. position dependent BCs). The necessarily value to indicate an
// outflow BC differs between advection and the NodalProj, so we pass it in
void incflo::prob_set_BC_MF (Orientation ori, Box const& bx,
                             Array4<int> const& mask, int lev,
                             int outflow_val)
{
    if (1100 == m_probtype || 1101 == m_probtype || 1102 == m_probtype)
    {
        int direction = 0;
        if (1101 == m_probtype) {
            direction = 1;
        }
        else if (1102 == m_probtype) {
            direction = 2;
        }
        Box const& domain = geom[lev].Domain();
        int half_num_cells  = domain.length(direction) / 2;

        Orientation::Side side = ori.faceDir();
        if (side == Orientation::low) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (direction == 0) {
                    if (i <= half_num_cells) {
                        mask(i,j,k,0) = outflow_val; // outflow on bottom
                    } else {
                        mask(i,j,k,0) = BCType::ext_dir; // inflow on top
                    }
                }
                else if (direction == 1) {
                    if (j <= half_num_cells) {
                        mask(i,j,k,0) = outflow_val;
                    } else {
                        mask(i,j,k,0) = BCType::ext_dir;
                    }
                }
                else if (direction == 2) {
                    if (k <= half_num_cells) {
                        mask(i,j,k,0) = outflow_val;
                    } else {
                        mask(i,j,k,0) = BCType::ext_dir;
                    }
                }
            });
        } else {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (direction == 0) {
                    if (i > half_num_cells) {
                        mask(i,j,k,0) = outflow_val; // outflow on top
                    } else {
                        mask(i,j,k,0) = BCType::ext_dir;  // inflow on bottom
                    }
                }
                else if (direction == 1) {
                    if (j > half_num_cells) {
                        mask(i,j,k,0) = outflow_val;
                    } else {
                        mask(i,j,k,0) = BCType::ext_dir;
                    }
                }
                else if (direction == 2) {
                    if (k > half_num_cells) {
                        mask(i,j,k,0) = outflow_val;
                    } else {
                        mask(i,j,k,0) = BCType::ext_dir;
                    }
                }
            });
        }
    }
    else
    {
        Abort("incflo::prob_set_BC_MF: No masking function for probtype "
              +std::to_string(m_probtype));
    }
}

// FIXME - need one like this for diffusion, but then it really needs to be going to the PhysBCFunct
// Because we need the Dirichelet val, right?
// For MAC, the BC is on phi, so always 0 or 1, I think
void incflo::prob_set_MAC_robinBCs (Orientation ori, Box const& bx,
                                    Array4<Real> const& robin_a,
                                    Array4<Real> const& robin_b,
                                    Array4<Real> const& robin_f,
                                    int lev)
{
    // Robin BC:   a u + b du/dn = f  -- inflow,  Neumann   a=0, b=1, f=0
    //                                -- outflow, Dirichlet a=1, b=0, f=0

    if (1100 == m_probtype || 1101 == m_probtype || 1102 == m_probtype)
    {
        int direction = 0;
        if (1101 == m_probtype) {
            direction = 1;
        }
        else if (1102 == m_probtype) {
            direction = 2;
        }
        Box const& domain = geom[lev].Domain();
        int half_num_cells  = domain.length(direction) / 2;

        Orientation::Side side = ori.faceDir();
        if (side == Orientation::low) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                robin_f(i,j,k) = 0.;

                if (direction == 0) {
                    if (i <= half_num_cells) {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0; // outflow on bottom
                    } else {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.; // inflow on top
                    }
                }
                else if (direction == 1) {
                    if (j <= half_num_cells) {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0;
                    } else {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.;
                    }
                }
                else if (direction == 2) {
                    if (k <= half_num_cells) {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0;
                    } else {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.;
                    }
                }
            });
        } else {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                robin_f(i,j,k,0) = 0.;

                if (direction == 0) {
                    if (i > half_num_cells) {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0; // outflow on top
                    } else {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.;  // inflow on bottom
                    }
                }
                else if (direction == 1) {
                    if (j > half_num_cells) {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0;
                    } else {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.;
                    }
                }
                else if (direction == 2) {
                    if (k > half_num_cells) {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0;
                    } else {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.;
                    }
                }
            });
        }
    }
    else
    {
        Abort("incflo::prob_set_MAC_robinBCs: No masking function for probtype "
              +std::to_string(m_probtype));
    }
}
