#include <incflo.H>

using namespace amrex;

// Specify mixed BCs for both the nodal projection and advection. The necessary value
// to indicate an in/outflow BC differs between advection and the NodalProj, so we pass it.
// Note that the advection BCs are on the velocity and scalars, whereas for the Nodal
// Projection we're creating an overset mask that operates on the solve for phi (~pressure)
void incflo::prob_set_BC_MF (Orientation const& ori, Box const& bx,
                             Array4<int> const& mask, int lev,
                             int inflow_val, int outflow_val,
                             std::string const& field)
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

        // for this problem, bcs are same for all fields, only ncomp varies
        int ncomp =  field == "velocity" ? AMREX_SPACEDIM : 1;

        Orientation::Side side = ori.faceDir();
        if (side == Orientation::low) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                for ( int n = 0; n < ncomp; n++ ){
                    if (direction == 0) {
                        if (i <= half_num_cells) {
                            mask(i,j,k,n) = outflow_val; // outflow on bottom
                        } else {
                            mask(i,j,k,n) = inflow_val; // inflow on top
                        }
                    }
                    else if (direction == 1) {
                        if (j <= half_num_cells) {
                            mask(i,j,k,n) = outflow_val;
                        } else {
                            mask(i,j,k,n) = inflow_val;
                        }
                    }
#if (AMREX_SPACEDIM == 3)
                    else if (direction == 2) {
                        if (k <= half_num_cells) {
                            mask(i,j,k,n) = outflow_val;
                        } else {
                            mask(i,j,k,n) = inflow_val;
                        }
                    }
#endif
                }
            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                for ( int n = 0; n < ncomp; n++ ){
                    if (direction == 0) {
                        if (i > half_num_cells) {
                            mask(i,j,k,n) = outflow_val; // outflow on top
                        } else {
                            mask(i,j,k,n) = inflow_val;  // inflow on bottom
                        }
                    }
                    else if (direction == 1) {
                        if (j > half_num_cells) {
                            mask(i,j,k,n) = outflow_val;
                        } else {
                            mask(i,j,k,n) = inflow_val;
                        }
                    }
#if (AMREX_SPACEDIM == 3)
                    else if (direction == 2) {
                        if (k > half_num_cells) {
                            mask(i,j,k,n) = outflow_val;
                        } else {
                            mask(i,j,k,n) = inflow_val;
                        }
                    }
#endif
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

void incflo::prob_set_MAC_robinBCs (Orientation const& ori, Box const& bx,
                                    Array4<Real> const& robin_a,
                                    Array4<Real> const& robin_b,
                                    Array4<Real> const& robin_f,
                                    int lev)
{
    // For the MAC Projection, the BC is on phi (~pressure), so always homogeneous
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
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                robin_f(i,j,k) = 0.;

                if (direction == 0) {
                    if (i <= half_num_cells) {
                        robin_a(i,j,k) = 1.;
                        robin_b(i,j,k) = 0; // outflow on bottom
                    } else {
                        robin_a(i,j,k) = 0.;
                        robin_b(i,j,k) = 1.; // inflow on top
                    }
                }
                else if (direction == 1) {
                    if (j <= half_num_cells) {
                        robin_a(i,j,k) = 1.;
                        robin_b(i,j,k) = 0;
                    } else {
                        robin_a(i,j,k) = 0.;
                        robin_b(i,j,k) = 1.;
                    }
                }
#if (AMREX_SPACEDIM == 3)
                else if (direction == 2) {
                    if (k <= half_num_cells) {
                        robin_a(i,j,k) = 1.;
                        robin_b(i,j,k) = 0;
                    } else {
                        robin_a(i,j,k) = 0.;
                        robin_b(i,j,k) = 1.;
                    }
                }
#endif
            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                robin_f(i,j,k) = 0.;

                if (direction == 0) {
                    if (i > half_num_cells) {
                        robin_a(i,j,k) = 1.;
                        robin_b(i,j,k) = 0; // outflow on top
                    } else {
                        robin_a(i,j,k) = 0.;
                        robin_b(i,j,k) = 1.;  // inflow on bottom
                    }
                }
                else if (direction == 1) {
                    if (j > half_num_cells) {
                        robin_a(i,j,k) = 1.;
                        robin_b(i,j,k) = 0;
                    } else {
                        robin_a(i,j,k) = 0.;
                        robin_b(i,j,k) = 1.;
                    }
                }
#if (AMREX_SPACEDIM == 3)
                else if (direction == 2) {
                    if (k > half_num_cells) {
                        robin_a(i,j,k) = 1.;
                        robin_b(i,j,k) = 0;
                    } else {
                        robin_a(i,j,k) = 0.;
                        robin_b(i,j,k) = 1.;
                    }
                }
#endif
            });
        }
    }
    else
    {
        Abort("incflo::prob_set_MAC_robinBCs: No masking function for probtype "
              +std::to_string(m_probtype));
    }
}

void incflo::prob_set_diffusion_robinBCs (Orientation const& ori, Box const& bx,
                                          Array4<Real> const& robin_a,
                                          Array4<Real> const& robin_b,
                                          Array4<Real> const& robin_f,
                                          Array4<Real const> const& bcval,
                                          int lev)
{
    // For diffusion, we also pass in the dirichlet bc (bcval)
    // Robin BC:   a u + b du/dn = f  -- inflow,  Dirichlet a=1, b=0, f=bcval
    //                                -- outflow, Neumann   a=0, b=1, f=0

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
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (direction == 0) {
                    if (i <= half_num_cells) {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.; // outflow on bottom
                        robin_f(i,j,k,0) = 0.;
                    } else {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0.; // inflow on top
                        robin_f(i,j,k,0) = bcval(i,j,k,0);
                    }
                }
                else if (direction == 1) {
                    if (j <= half_num_cells) {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.;
                        robin_f(i,j,k,0) = 0.;
                    } else {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0.;
                        robin_f(i,j,k,0) = bcval(i,j,k,0);
                    }
                }
#if (AMREX_SPACEDIM == 3)
                else if (direction == 2) {
                    if (k <= half_num_cells) {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.;
                        robin_f(i,j,k,0) = 0.;
                    } else {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0.;
                        robin_f(i,j,k,0) = bcval(i,j,k,0);
                    }
                }
#endif
            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                robin_f(i,j,k,0) = 0.;

                if (direction == 0) {
                    if (i > half_num_cells) {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.; // outflow on top
                        robin_f(i,j,k,0) = 0.;
                    } else {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0.;  // inflow on bottom
                        robin_f(i,j,k,0) = bcval(i,j,k,0);
                    }
                }
                else if (direction == 1) {
                    if (j > half_num_cells) {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.;
                        robin_f(i,j,k,0) = 0.;
                    } else {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0.;
                        robin_f(i,j,k,0) = bcval(i,j,k,0);
                    }
                }
#if (AMREX_SPACEDIM == 3)
                else if (direction == 2) {
                    if (k > half_num_cells) {
                        robin_a(i,j,k,0) = 0.;
                        robin_b(i,j,k,0) = 1.;
                        robin_f(i,j,k,0) = 0.;
                    } else {
                        robin_a(i,j,k,0) = 1.;
                        robin_b(i,j,k,0) = 0.;
                        robin_f(i,j,k,0) = bcval(i,j,k,0);
                    }
                }
#endif
            });
        }
    }
    else
    {
        Abort("incflo::prob_set_diffusion_robinBCs: No masking function for probtype "
              +std::to_string(m_probtype));
    }
}

void incflo::prob_set_inflow_velocity (int /*grid_id*/, Orientation ori, Box const& bx,
                                       Array4<Real> const& vel, int lev, Real /*time*/)
{
    if (6 == m_probtype)
    {
        AMREX_D_TERM(Real u = m_ic_u;,
                     Real v = m_ic_v;,
                     Real w = m_ic_w;);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            AMREX_D_TERM(vel(i,j,k,0) = u;,
                         vel(i,j,k,1) = v;,
                         vel(i,j,k,2) = w;);
        });
    }
    else if (31 == m_probtype)
    {
        Real dyinv = 1.0 / Geom(lev).Domain().length(1);
        Real u = m_ic_u;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real y = (j+0.5)*dyinv;
            vel(i,j,k,0) = 6. * u * y * (1.-y);
        });
    }
#if (AMREX_SPACEDIM == 3)
    else if (311 == m_probtype)
    {
        Real dzinv = 1.0 / Geom(lev).Domain().length(2);
        Real u = m_ic_u;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = (k+0.5)*dzinv;
            vel(i,j,k,0) = 6. * u * z * (1.-z);
        });
    }
    else if (41 == m_probtype)
    {
        Real dzinv = 1.0 / Geom(lev).Domain().length(2);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = (k+0.5)*dzinv;
            vel(i,j,k,0) = 0.5*z;
        });
    }
    else if (32 == m_probtype)
    {
        Real dzinv = 1.0 / Geom(lev).Domain().length(2);
        Real v = m_ic_v;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = (k+0.5)*dzinv;
            vel(i,j,k,1) = 6. * v * z * (1.-z);
        });
    }
#endif
    else if (322 == m_probtype)
    {
        Real dxinv = 1.0 / Geom(lev).Domain().length(0);
        Real v = m_ic_v;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = (i+0.5)*dxinv;
            vel(i,j,k,1) = 6. * v * x * (1.-x);
        });
    }
    else if (33 == m_probtype)
    {
        Real dxinv = 1.0 / Geom(lev).Domain().length(0);
        Real w = m_ic_w;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = (i+0.5)*dxinv;
            vel(i,j,k,2) = 6. * w * x * (1.-x);
        });
    }
    else if (333 == m_probtype)
    {
        Real dyinv = 1.0 / Geom(lev).Domain().length(1);
        Real w = m_ic_w;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real y = (j+0.5)*dyinv;
            vel(i,j,k,2) = 6. * w * y * (1.-y);
        });
    }
    else
    {
        const int  dir = ori.coordDir();
        const Real bcv = m_bc_velocity[ori][dir];
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            vel(i,j,k,dir) = bcv;
        });
    };
}
