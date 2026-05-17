#include <incflo.H>

using namespace amrex;

// Compute temperature forcing terms.
void incflo::compute_tem_forces (Real /*time*/, Vector<MultiFab*> const& tem_forces)
{
    if (m_use_temperature) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (MFIter mfi(*tem_forces[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& tem_f = tem_forces[lev]->array(mfi);

                ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // Just H_T here, as in
                    //   rho*cp ( dT/dt + U dot grad T) = div mu_T grad T + H_T
                    // For now we don't have any external forces on temperature
                    tem_f(i,j,k) = 0.0;
                });
            }
        }
    }
}

void incflo::compute_tra_forces (Vector<MultiFab*> const& tra_forces,
                                 Vector<MultiFab const*> const& density)
{
    if (m_advect_tracer) {

        auto const* iconserv = get_tracer_iconserv_device_ptr();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (MFIter mfi(*tra_forces[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real>       const& tra_f = tra_forces[lev]->array(mfi);
                Array4<Real const> const& rho   =    density[lev]->const_array(mfi);

                ParallelFor(bx, m_ntrac,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    // For now we don't have any external forces on the scalars
                    tra_f(i,j,k,n) = 0.0;

                    if (iconserv[n]){
                        // Return the force term for the update of (rho s), NOT just s.
                        tra_f(i,j,k,n) *= rho(i,j,k);
                    }
                });
            }
        }
    }
}

void incflo::compute_vel_forces (Vector<MultiFab*> const& vel_forces,
                                 Vector<MultiFab const*> const& velocity,
                                 Vector<MultiFab const*> const& density,
                                 Vector<MultiFab const*> const& tracer_old,
                                 Vector<MultiFab const*> const& tracer_new,
                                 bool include_pressure_gradient)
{
    for (int lev = 0; lev <= finest_level; ++lev)
       compute_vel_forces_on_level (lev, *vel_forces[lev], *velocity[lev], *density[lev],
                                         *tracer_old[lev], *tracer_new[lev], include_pressure_gradient);
}

void incflo::compute_vel_forces_on_level (int lev,
                                                MultiFab& vel_forces,
                                          const MultiFab& /*velocity*/,
                                          const MultiFab& density,
                                          const MultiFab& tracer_old,
                                          const MultiFab& tracer_new,
                                          bool include_pressure_gradient)
{
    GpuArray<Real,3> l_gravity{m_gravity[0],m_gravity[1],m_gravity[2]};
    GpuArray<Real,3> l_gp0{m_gp0[0], m_gp0[1], m_gp0[2]};

    auto const dx = geom[lev].CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
            Box const& bx = mfi.tilebox();
            Array4<Real>       const& vel_f =  vel_forces.array(mfi);
            Array4<Real const> const&   rho =     density.const_array(mfi);
            Array4<Real const> const& gradp = m_leveldata[lev]->gp.const_array(mfi);

            if (m_use_boussinesq) {
                // This uses a Boussinesq approximation where the buoyancy depends on
                //      first tracer rather than density
                Array4<Real const> const& tra_o = tracer_old.const_array(mfi);
                Array4<Real const> const& tra_n = tracer_new.const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    int n = 0; // Potential temperature

                    Real rhoinv = Real(1.0)/rho(i,j,k);
                    Real ft = Real(0.5) * (tra_o(i,j,k,n) + tra_n(i,j,k,n));

                    if (include_pressure_gradient)
                    {
                        AMREX_D_TERM(vel_f(i,j,k,0) = -gradp(i,j,k,0)*rhoinv + l_gravity[0] * ft;,
                                     vel_f(i,j,k,1) = -gradp(i,j,k,1)*rhoinv + l_gravity[1] * ft;,
                                     vel_f(i,j,k,2) = -gradp(i,j,k,2)*rhoinv + l_gravity[2] * ft;);
                    } else {
                        AMREX_D_TERM(vel_f(i,j,k,0) =                          l_gravity[0] * ft;,
                                     vel_f(i,j,k,1) =                          l_gravity[1] * ft;,
                                     vel_f(i,j,k,2) =                          l_gravity[2] * ft;);
                    }
                });

            } else if (m_probtype == 16) {
                Real Re = Real(1)/m_mu;  // Note this assumes you are running exactly the problem set up, with U = 1 and L = 1 and rho = 1.
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhoinv = Real(1)/rho(i,j,k);

                    if (include_pressure_gradient)
                    {
                        AMREX_D_TERM(vel_f(i,j,k,0) = -(gradp(i,j,k,0)+l_gp0[0])*rhoinv + l_gravity[0];,
                                     vel_f(i,j,k,1) = -(gradp(i,j,k,1)+l_gp0[1])*rhoinv + l_gravity[1];,
                                     vel_f(i,j,k,2) = -(gradp(i,j,k,2)+l_gp0[2])*rhoinv + l_gravity[2];);
                    } else {
                        AMREX_D_TERM(vel_f(i,j,k,0) = -(               l_gp0[0])*rhoinv + l_gravity[0];,
                                     vel_f(i,j,k,1) = -(               l_gp0[1])*rhoinv + l_gravity[1];,
                                     vel_f(i,j,k,2) = -(               l_gp0[2])*rhoinv + l_gravity[2];);
                    }

                    Real x = (i+Real(0.5)) * dx[0];
                    Real y = (j+Real(0.5)) * dx[1];

                    Real f     = x*x*x*x - Real(2)*x*x*x + x*x;
                    Real g     = y*y*y*y - y*y;
                    Real capF  = Real(0.2) * x*x*x*x*x   - Real(0.5) * x*x*x*x   + (Real(1)/Real(3)) * x*x*x;
                    Real capF1 = Real(-4)  * x*x*x*x*x*x + Real(12)  * x*x*x*x*x -  Real(14)         * x*x*x*x + Real(8) * x*x*x - Real(2) * x*x;
                    Real capF2 = Real(0.5) * f * f;
                    Real capG1 = Real(-24) * y*y*y*y*y + Real(8) * y*y*y - Real(4) * y;

                    Real  fp   = Real(4) * x*x*x - Real(6)*x*x + Real(2)*x;
                    //Real  fpp  = Real(12) * x*x - Real(12)*x + Real(2);
                    Real  fppp = Real(24) * x - Real(12);

                    Real  gp   =  Real(4) * y*y*y - Real(2)*y;
                    Real  gpp  = Real(12) * y*y - Real(2);

                    vel_f(i,j,k,1) += Real(8) / Re * (Real(24) * capF + Real(2) * fp * gpp + fppp * g) + Real(64) * (capF2 * capG1 - g * gp * capF1);
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhoinv = Real(1)/rho(i,j,k);

                    if (include_pressure_gradient)
                    {
                        AMREX_D_TERM(vel_f(i,j,k,0) = -(gradp(i,j,k,0)+l_gp0[0])*rhoinv + l_gravity[0];,
                                     vel_f(i,j,k,1) = -(gradp(i,j,k,1)+l_gp0[1])*rhoinv + l_gravity[1];,
                                     vel_f(i,j,k,2) = -(gradp(i,j,k,2)+l_gp0[2])*rhoinv + l_gravity[2];);
                    } else {
                        AMREX_D_TERM(vel_f(i,j,k,0) = -(               l_gp0[0])*rhoinv + l_gravity[0];,
                                     vel_f(i,j,k,1) = -(               l_gp0[1])*rhoinv + l_gravity[1];,
                                     vel_f(i,j,k,2) = -(               l_gp0[2])*rhoinv + l_gravity[2];);
                    }
                });
            }
    }
}
