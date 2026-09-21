#include <incflo.H>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace amrex;

//
// Compute new dt by using the formula derived in
// "A Boundary Condition Capturing Method for Multiphase Incompressible Flow"
// by Kang et al. (JCP).
//
//  dt/2 * ( C+V + sqrt( (C+V)**2 + 4Fx/dx + 4Fy/dy + 4Fz/dz )
//
// where
//
// C = max(|U|)/dx + max(|V|)/dy + max(|W|)/dz    --> Convection
//
// V = 2 * max(eta/rho) * (1/dx^2 + 1/dy^2 +1/dz^2) --> Diffusion
//
// Fx, Fy, Fz = net acceleration due to external forces
//
// WARNING: We use a slightly modified version of C in the implementation below
//
void incflo::ComputeDt (int initialization, bool explicit_diffusion, double cur_time)
{
    BL_PROFILE("incflo::ComputeDt");

    // Store the past two dt
    m_prev_prev_dt = m_prev_dt;
    m_prev_dt = m_dt;

    Real conv_cfl = Real(0.0);
    Real diff_cfl = Real(0.0);
    Real forc_cfl = Real(0.0);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto const dxinv = geom[lev].InvCellSizeArray();
        MultiFab const& vel   = m_leveldata[lev]->velocity;
        MultiFab const& rho   = m_leveldata[lev]->density;
        MultiFab const& tra   = m_leveldata[lev]->tracer;
        MultiFab const& tra_o = m_leveldata[lev]->tracer_o;

        Real conv_lev = Real(0.0);
        Real diff_lev = Real(0.0);
        Real forc_lev = Real(0.0);

       // Make a temporary here to hold vel_forces
       MultiFab vel_forces(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);

       compute_vel_forces_on_level (lev, vel_forces, vel, rho, tra_o, tra);

       // Explicit-diffusion bound: the largest diffusivity that is applied
       // explicitly, divided by rho.  We must use the strain-rate dependent
       // viscosity here (not just the constant m_mu) and we must cover the tracer
       // and temperature diffusivities as well, since those are updated explicitly
       // with the same m_diff_type switch.  Covered cells are set to zero.
       MultiFab nu;
       if (explicit_diffusion) {
           nu.define(grids[lev], dmap[lev], 1, 0, MFInfo(), Factory(lev));
           compute_viscosity_at_level(lev, &nu, &m_leveldata[lev]->density,
                                      &m_leveldata[lev]->velocity, geom[lev],
                                      m_cur_time, 0);
           Real mu_s_max = Real(0.0);
           if (m_advect_tracer) {
               for (int n = 0; n < m_ntrac; ++n) {
                   mu_s_max = amrex::max(mu_s_max, m_mu_s[n]);
               }
           }
           Real mu_T_eff = m_use_temperature ? m_mu_T / m_cp : Real(0.0);
#ifdef AMREX_USE_EB
           auto const& dt_flags = EBFactory(lev).getMultiEBCellFlagFab();
#endif
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
           for (MFIter mfi(nu,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
               Box const& bx = mfi.tilebox();
               Array4<Real> const& nu_a = nu.array(mfi);
               Array4<Real const> const& r = rho.const_array(mfi);
#ifdef AMREX_USE_EB
               Array4<EBCellFlag const> const& f = dt_flags.const_array(mfi);
#endif
               ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
#ifdef AMREX_USE_EB
                   if (f(i,j,k).isCovered()) { nu_a(i,j,k) = Real(0.0); return; }
#endif
                   Real rinv = Real(1.0)/r(i,j,k);
                   // Velocity and temperature diffuse with eta/rho and mu_T/(rho cp);
                   // a non-conservative tracer diffuses with mu_s itself.
                   nu_a(i,j,k) = amrex::max(amrex::max(nu_a(i,j,k), mu_T_eff)*rinv,
                                            mu_s_max*amrex::max(Real(1.0), rinv));
               });
           }
       }

#ifdef AMREX_USE_EB
        if (!vel.isAllRegular()) {
            auto const& flag = EBFactory(lev).getMultiEBCellFlagFab();
            conv_lev = amrex::ReduceMax(vel, flag, 0,
                       [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                                  Array4<Real const> const& v,
                                                  Array4<EBCellFlag const> const& f) -> Real
                       {
                           Real mx = -1.0;
                           amrex::Loop(b, [=,&mx] (int i, int j, int k) noexcept
                           {
                               if (!f(i,j,k).isCovered()) {
                                   mx = amrex::max(AMREX_D_DECL(amrex::Math::abs(v(i,j,k,0))*dxinv[0],
                                                                amrex::Math::abs(v(i,j,k,1))*dxinv[1],
                                                                amrex::Math::abs(v(i,j,k,2))*dxinv[2]), mx);
                               }
                           });
                           return mx;
                       });
            if (explicit_diffusion) {
                diff_lev = nu.max(0, 0, true);
            }

            // Forcing term -- old way of computing
            // const auto dxinv_finest = Geom(finest_level).InvCellSizeArray();
            // forc_lev = std::abs(m_gravity[0] - std::abs(m_gp0[0])) * dxinv_finest[0]
            //          + std::abs(m_gravity[1] - std::abs(m_gp0[1])) * dxinv_finest[1]
            //          + std::abs(m_gravity[2] - std::abs(m_gp0[2])) * dxinv_finest[2];

            // Forcing term -- new way of computing using "actual" forcing term
            forc_lev = amrex::ReduceMax(vel_forces, flag, 0,
                  [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                             Array4<Real const> const& vf,
                                             Array4<EBCellFlag const> const& f) -> Real
                  {
                      Real mx = Real(-1.0);
                      amrex::Loop(b, [=,&mx] (int i, int j, int k) noexcept
                      {
                          if (!f(i,j,k).isCovered()) {
                              mx = amrex::max(AMREX_D_DECL(amrex::Math::abs(vf(i,j,k,0))*dxinv[0],
                                                           amrex::Math::abs(vf(i,j,k,1))*dxinv[1],
                                                           amrex::Math::abs(vf(i,j,k,2))*dxinv[2]), mx);
                          }
                      });
                      return mx;
                  });
        } else
#endif
        {
            conv_lev = amrex::ReduceMax(vel, 0,
                       [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                                  Array4<Real const> const& v) -> Real
                       {
                           Real mx = Real(-1.0);
                           amrex::Loop(b, [=,&mx] (int i, int j, int k) noexcept
                           {
                               mx = amrex::max(AMREX_D_DECL(amrex::Math::abs(v(i,j,k,0))*dxinv[0],
                                                            amrex::Math::abs(v(i,j,k,1))*dxinv[1],
                                                            amrex::Math::abs(v(i,j,k,2))*dxinv[2]), mx);
                           });
                           return mx;
                       });

            if (explicit_diffusion) {
                diff_lev = nu.max(0, 0, true);
            }

            // Forcing term -- old way of computing
            // const auto dxinv_finest = Geom(finest_level).InvCellSizeArray();
            // forc_lev = std::abs(m_gravity[0] - std::abs(m_gp0[0])) * dxinv_finest[0]
            //          + std::abs(m_gravity[1] - std::abs(m_gp0[1])) * dxinv_finest[1]
            //          + std::abs(m_gravity[2] - std::abs(m_gp0[2])) * dxinv_finest[2];

            // Forcing term -- new way of computing using "actual" forcing term
            forc_lev = amrex::ReduceMax(vel_forces, 0,
                  [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                             Array4<Real const> const& vf) -> Real
                  {
                      Real mx = Real(-1.0);
                      amrex::Loop(b, [=,&mx] (int i, int j, int k) noexcept
                      {
                          mx = amrex::max(AMREX_D_DECL(amrex::Math::abs(vf(i,j,k,0))*dxinv[0],
                                                       amrex::Math::abs(vf(i,j,k,1))*dxinv[1],
                                                       amrex::Math::abs(vf(i,j,k,2))*dxinv[2]), mx);
                      });
                      return mx;
                  });
        }

        forc_cfl = std::max(forc_cfl, forc_lev);
        conv_cfl = std::max(conv_cfl, conv_lev);

#if (AMREX_SPACEDIM == 2)
        Real dxinv_norm = dxinv[0]*dxinv[0]+dxinv[1]*dxinv[1];
#else
        Real dxinv_norm = dxinv[0]*dxinv[0]+dxinv[1]*dxinv[1]+dxinv[2]*dxinv[2];
#endif

        diff_cfl = std::max(diff_cfl, diff_lev*Real(2.0)*dxinv_norm);
    }

    Real cd_cfl;
    if (explicit_diffusion) {
        ParallelAllReduce::Max<Real>({conv_cfl,diff_cfl},
                                     ParallelContext::CommunicatorSub());
        cd_cfl = conv_cfl + diff_cfl;
    } else {
        ParallelAllReduce::Max<Real>(conv_cfl,
                                     ParallelContext::CommunicatorSub());
        cd_cfl = conv_cfl;
    }

    ParallelAllReduce::Max<Real>(forc_cfl,
                                 ParallelContext::CommunicatorSub());

    // Combined CFL conditioner
    Real comb_cfl = cd_cfl + std::sqrt(cd_cfl*cd_cfl + Real(4.0) * forc_cfl);

    // Update dt
    Real dt_new;
    if (comb_cfl > 0.)
    {
        dt_new = Real(2.0) * m_cfl / comb_cfl;

    } else {

        // This is totally random but just a way to set a timestep
        // when the initial velocity is zero and the forcing term
        // is not a body force
        auto const dx    = geom[finest_level].CellSizeArray();
        dt_new = std::min(dx[0],dx[1]);
#if (AMREX_SPACEDIM == 3)
        dt_new = std::min(dt_new,dx[2]);
#endif
    }

    // Optionally reduce CFL for initial step
    if(initialization)
    {
        dt_new *= m_init_shrink;
    }

    // Protect against very small comb_cfl
    // This may happen, for example, when the initial velocity field
    // is zero for an inviscid flow with no external forcing
    Real eps = std::numeric_limits<Real>::epsilon();
    if(! initialization && comb_cfl <= eps)
    {
        dt_new = Real(0.5) * m_dt;
    }

    // Don't let the timestep grow by more than m_dt_change_max per step
    // unless the previous time step was unduly shrunk to match m_plot_per_exact
    Real allowed_change_factor = m_dt_change_max;
    if( (m_dt > Real(0.0)) && !(m_plot_per_exact > 0 && m_last_plt == m_nstep && m_nstep > 0) )
    {
        dt_new = amrex::min(dt_new, allowed_change_factor * m_prev_dt);
    }
    else if ( (m_dt > Real(0.0)) && (m_plot_per_exact > 0 && m_last_plt == m_nstep && m_nstep > 0) )
    {
        dt_new = amrex::min( dt_new, allowed_change_factor * amrex::max(m_prev_dt, m_prev_prev_dt) );
    }

    // Do not overshoot specified plot times. Use double-precision cur_time
    // so single-precision builds do not clip because of accumulated time drift.
    if (m_plot_per_exact > Real(0.0))
    {
        double const plot_per_exact = m_plot_per_exact;
        double const dt_new_d = dt_new;
        double const eps_d = eps;
        if (std::trunc((cur_time + dt_new_d + eps_d) / plot_per_exact) >
            std::trunc((cur_time + eps_d) / plot_per_exact))
        {
            dt_new = static_cast<Real>(
                std::trunc((cur_time + dt_new_d) / plot_per_exact) * plot_per_exact - cur_time);
        }
    }

    // Do not overshoot the final time if not running to steady state.
    if (!m_steady_state && m_stop_time > Real(0.0))
    {
        double const stop_time = m_stop_time;
        double const dt_new_d = dt_new;
        if (cur_time + dt_new_d > stop_time)
        {
            dt_new = static_cast<Real>(stop_time - cur_time);
        }
    }

    // Make sure the timestep is not set to zero after a m_plot_per_exact stop
    if (dt_new < eps)
    {
        dt_new = Real(0.5) * m_dt;
    }

    // If using fixed time step, check CFL condition and give warning if not satisfied
    if (m_fixed_dt > Real(0.0))
    {
        if(dt_new < m_fixed_dt)
        {
            amrex::Print() << "WARNING: fixed_dt does not satisfy CFL condition: \n"
                           << "max dt by CFL     : " << dt_new << "\n"
                           << "fixed dt specified: " << m_fixed_dt << "\n";
        }
        m_dt = m_fixed_dt;
    }
    else
    {
        m_dt = dt_new;
    }
}
