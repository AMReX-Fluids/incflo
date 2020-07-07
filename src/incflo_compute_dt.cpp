#include <incflo.H>

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
void incflo::ComputeDt (int initialization, bool explicit_diffusion)
{
    BL_PROFILE("incflo::ComputeDt");

    // Store the past two dt
    m_prev_prev_dt = m_prev_dt;
    m_prev_dt = m_dt;

    Real conv_cfl = 0.0;
    Real diff_cfl = 0.0;
    Real forc_cfl = 0.0;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto const dxinv = geom[lev].InvCellSizeArray();
        MultiFab const& vel   = m_leveldata[lev]->velocity;
        MultiFab const& rho   = m_leveldata[lev]->density;
        MultiFab const& tra   = m_leveldata[lev]->tracer;
        MultiFab const& tra_o = m_leveldata[lev]->tracer_o;

        Real conv_lev = 0.0;
        Real diff_lev = 0.0;
        Real forc_lev = 0.0;

       // Make a temporary here to hold vel_forces
       MultiFab vel_forces(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);

       compute_vel_forces_on_level (lev, vel_forces, vel, rho, tra_o, tra);

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
                diff_lev = amrex::ReduceMax(rho, flag, 0,
                           [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                                      Array4<Real const> const& r,
                                                      Array4<EBCellFlag const> const& f) -> Real
                          {
                              Real mx = -1.0;
                              amrex::Loop(b, [=,&mx] (int i, int j, int k) noexcept
                              {
                                  if (!f(i,j,k).isCovered()) {
                                      Real rho_inv = 1.0/r(i,j,k);
                                      mx = amrex::max(rho_inv, mx);
                                  }
                              });
                              return mx;
                          });
                diff_lev *= m_mu;
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
                      Real mx = -1.0;
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
                           Real mx = -1.0;
                           amrex::Loop(b, [=,&mx] (int i, int j, int k) noexcept
                           {
                               mx = amrex::max(AMREX_D_DECL(amrex::Math::abs(v(i,j,k,0))*dxinv[0],
                                                            amrex::Math::abs(v(i,j,k,1))*dxinv[1],
                                                            amrex::Math::abs(v(i,j,k,2))*dxinv[2]), mx);
                           });
                           return mx;
                       });

            if (explicit_diffusion) {
                diff_lev = amrex::ReduceMax(rho, 0,
                           [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                                      Array4<Real const> const& r) -> Real
                           {
                               Real mx = -1.0;
                               amrex::Loop(b, [=,&mx] (int i, int j, int k) noexcept
                               {
                                   Real rho_inv = 1.0/r(i,j,k);
                                   mx = amrex::max(rho_inv, mx);
                               });
                               return mx;
                           });
                diff_lev *= m_mu;
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
                      Real mx = -1.0;
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
        diff_cfl = std::max(diff_cfl, diff_lev*2.0_rt*(dxinv[0]*dxinv[0]+dxinv[1]*dxinv[1]+
                                                       dxinv[2]*dxinv[2]));
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
    Real comb_cfl = cd_cfl + std::sqrt(cd_cfl*cd_cfl + 4.0 * forc_cfl);

    // Update dt
    Real dt_new = 2.0 * m_cfl / comb_cfl;

    // Optionally reduce CFL for initial step
    if(initialization)
    {
        dt_new *= m_init_shrink;
    }

    // Protect against very small comb_cfl
    // This may happen, for example, when the initial velocity field
    // is zero for an inviscid flow with no external forcing
    Real eps = std::numeric_limits<Real>::epsilon();
    if(comb_cfl <= eps)
    {
        dt_new = 0.5 * m_dt;
    }

    // Don't let the timestep grow by more than 10% per step 
    // unless the previous time step was unduly shrunk to match m_plot_per_exact
    Real allowed_change_factor = 1.1;
    if( (m_dt > 0.0) && !(m_plot_per_exact > 0 && m_last_plt == m_nstep && m_nstep > 0) )
    {
        dt_new = amrex::min(dt_new, allowed_change_factor * m_prev_dt);
    } 
    else if ( (m_dt > 0.0) && (m_plot_per_exact > 0 && m_last_plt == m_nstep && m_nstep > 0) )
    {
        dt_new = amrex::min( dt_new, allowed_change_factor * amrex::max(m_prev_dt, m_prev_prev_dt) );
    }
    
    // Don't overshoot specified plot times
    if(m_plot_per_exact > 0.0 && 
            (std::trunc((m_cur_time + dt_new + eps) / m_plot_per_exact) > std::trunc((m_cur_time + eps) / m_plot_per_exact)))
    {
        dt_new = std::trunc((m_cur_time + dt_new) / m_plot_per_exact) * m_plot_per_exact - m_cur_time;
    }

    // Don't overshoot the final time if not running to steady state
    if(!m_steady_state && m_stop_time > 0.0)
    {
        if(m_cur_time + dt_new > m_stop_time)
        {
            dt_new = m_stop_time - m_cur_time;
        }
    }

    // Make sure the timestep is not set to zero after a m_plot_per_exact stop
    if(dt_new < eps)
    {
        dt_new = 0.5 * m_dt;
    }

    // If using fixed time step, check CFL condition and give warning if not satisfied
    if(m_fixed_dt > 0.0)
    {
	if(dt_new < m_fixed_dt)
	{
		amrex::Print() << "WARNING: fixed_dt does not satisfy CFL condition: \n"
					   << "max dt by CFL     : " << dt_new << "\n"
					   << "fixed dt specified: " << m_fixed_dt << std::endl;
	}
	m_dt = m_fixed_dt;
    }
    else
    {
	m_dt = dt_new;
    }
}
