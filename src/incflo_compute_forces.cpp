#include <incflo.H>

using namespace amrex;

void incflo::compute_tra_forces (Vector<MultiFab*> const& tra_forces,
                                 Vector<MultiFab const*> const& density)
{
    // NOTE: this routine must return the force term for the update of (rho s), NOT just s.
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
                Real Re = 1./m_mu;  // Note this assumes you are running exactly the problem set up, with U = 1 and L = 1 and rho = 1.
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhoinv = Real(1.0)/rho(i,j,k);

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

                    Real x = (i+0.5) * dx[0];
                    Real y = (j+0.5) * dx[1];

                    Real f     = x*x*x*x - 2.*x*x*x + x*x;
                    Real g     = y*y*y*y - y*y;
                    Real capF  =  0.2 * x*x*x*x*x   -  0.5 * x*x*x*x   + (1./3.)* x*x*x;
                    Real capF1 = -4.0 * x*x*x*x*x*x + 12.0 * x*x*x*x*x - 14.    * x*x*x*x + 8.0 * x*x*x - 2.0 * x*x;
                    Real capF2 = 0.5 * f * f;
                    Real capG1 = -24.0 * y*y*y*y*y + 8.0 * y*y*y - 4.0 * y;

                    Real  fp   = 4.0 * x*x*x - 6.0*x*x + 2.0*x;
                    //Real  fpp  = 12.0 * x*x - 12.0*x + 2.0;
                    Real  fppp = 24.0 * x - 12.0;

                    Real  gp   =  4.0 * y*y*y - 2.0*y;
                    Real  gpp  = 12.0 * y*y - 2.0;

                    vel_f(i,j,k,1) += 8.0 / Re * (24.0 * capF + 2.0 * fp * gpp + fppp * g) + 64.0 * (capF2 * capG1 - g * gp * capF1);
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhoinv = Real(1.0)/rho(i,j,k);

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
