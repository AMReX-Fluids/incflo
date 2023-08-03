#include <incflo.H>

using namespace amrex;

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
        Real dyinv = Real(1.0) / Geom(lev).Domain().length(1);
        Real u = m_ic_u;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real y = Real(j+0.5)*dyinv;
            vel(i,j,k,0) = Real(6.0) * u * y * (Real(1.0)-y);
        });
    }
    else if (311 == m_probtype)
    {
        Real dzinv = Real(1.0) / Geom(lev).Domain().length(2);
        Real u = m_ic_u;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = Real(k+0.5)*dzinv;
            vel(i,j,k,0) = Real(6.0) * u * z * (Real(1.0)-z);
        });
    }
    else if (41 == m_probtype)
    {
        Real dzinv = Real(1.0) / Geom(lev).Domain().length(2);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = Real(k+0.5)*dzinv;
            vel(i,j,k,0) = Real(0.5)*z;
        });
    }
    else if (32 == m_probtype)
    {
        Real dzinv = Real(1.0) / Geom(lev).Domain().length(2);
        Real v = m_ic_v;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = Real(k+0.5)*dzinv;
            vel(i,j,k,1) = Real(6.0) * v * z * (Real(1.0)-z);
        });
    }
    else if (322 == m_probtype)
    {
        Real dxinv = Real(1.0) / Geom(lev).Domain().length(0);
        Real v = m_ic_v;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = Real(i+0.5)*dxinv;
            vel(i,j,k,1) = Real(6.0) * v * x * (Real(1.0)-x);
        });
    }
    else if (33 == m_probtype)
    {
        Real dxinv = Real(1.0) / Geom(lev).Domain().length(0);
        Real w = m_ic_w;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = Real(i+0.5)*dxinv;
            vel(i,j,k,2) = Real(6.0) * w * x * (Real(1.0)-x);
        });
    }
    else if (333 == m_probtype)
    {
        Real dyinv = Real(1.0) / Geom(lev).Domain().length(1);
        Real w = m_ic_w;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real y = Real(j+0.5)*dyinv;
            vel(i,j,k,2) = Real(6.0) * w * y * (Real(1.0)-y);
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
