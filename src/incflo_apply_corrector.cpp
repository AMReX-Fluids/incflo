#include <incflo.H>

using namespace amrex;
//
// Apply corrector:
//
//  Output variables from the predictor are labelled _pred
//
//  1. Use u = vel_pred to compute
//
//      if (!advect_momentum) then
//          conv_u  = - u grad u
//      else
//          conv_u  = - del dot (rho u u)
//      conv_r  = - div( u rho  )
//      conv_t  = - div( u trac )
//      eta     = viscosity
//      divtau  = div( eta ( (grad u) + (grad u)^T ) ) / rho
//
//      conv_u  = 0.5 (conv_u + conv_u_pred)
//      conv_r  = 0.5 (conv_r + conv_r_pred)
//      conv_t  = 0.5 (conv_t + conv_t_pred)
//      if (m_diff_type == DiffusionType::Explicit)
//         divtau  = divtau at new_time using (*) state
//      else
//         divtau  = 0.0
//      eta     = eta at new_time
//
//     rhs = u + dt * ( conv + divtau )
//
//  2. Add explicit forcing term i.e. gravity + lagged pressure gradient
//
//      rhs += dt * ( g - grad(p + p0) / rho )
//
//  3. A. If (m_diff_type == DiffusionType::Implicit)
//        solve implicit diffusion equation for u*
//
//     ( 1 - dt / rho * div ( eta grad ) ) u* = u^n + dt * conv_u
//                                                  + dt * ( g - grad(p + p0) / rho )
//
//     B. If (m_diff_type == DiffusionType::Crank-Nicolson)
//        solve semi-implicit diffusion equation for u*
//
//     ( 1 - (dt/2) / rho * div ( eta grad ) ) u* = u^n + dt * conv_u + (dt/2) / rho * div (eta_old grad) u^n
//                                                      + dt * ( g - grad(p + p0) / rho )
//
//  4. Apply projection
//
//     Add pressure gradient term back to u*:
//
//      if (advect_momentum) then
//          (rho^(n+1) u**) = (rho^(n+1) u*) + dt * grad p
//      else
//          u** = u* + dt * grad p / rho^nph
//
//     Solve Poisson equation for phi:
//
//     div( grad(phi) / rho ) = div( u** )
//
//     Update pressure:
//
//     p = phi / dt
//
//     Update velocity, now divergence free
//
//     vel = u** - dt * grad p / rho
//
void incflo::ApplyCorrector()
{
    BL_PROFILE("incflo::ApplyCorrector");

    // We use the new time value for things computed on the "*" state
    Real new_time = m_cur_time + m_dt;

    // *************************************************************************************
    // Allocate space for the MAC velocities
    // *************************************************************************************
    Vector<MultiFab> AMREX_D_DECL(u_mac(finest_level+1), v_mac(finest_level+1), w_mac(finest_level+1));
    int ngmac = nghost_mac();

    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(u_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(0)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev));,
                     v_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(1)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev));,
                     w_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(2)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev)););
        if (ngmac > 0) {
            AMREX_D_TERM(u_mac[lev].setBndry(0.0);,
                         v_mac[lev].setBndry(0.0);,
                         w_mac[lev].setBndry(0.0););
        }
    }

    // **********************************************************************************************
    // We only reach the corrector if advection_type == MOL which means we don't use the forces
    //    in constructing the advection term
    // **********************************************************************************************
    Vector<MultiFab> vel_forces, tra_forces;
    Vector<MultiFab> vel_eta, tra_eta;
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel_forces.emplace_back(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost_force(),
                                MFInfo(), Factory(lev));
        if (m_advect_tracer) {
            tra_forces.emplace_back(grids[lev], dmap[lev], m_ntrac, nghost_force(),
                                    MFInfo(), Factory(lev));
        }
        vel_eta.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(), Factory(lev));
        if (m_advect_tracer) {
            tra_eta.emplace_back(grids[lev], dmap[lev], m_ntrac, 1, MFInfo(), Factory(lev));
        }
    }

    // *************************************************************************************
    // Compute the MAC-projected velocities at all levels
    // *************************************************************************************
    bool include_pressure_gradient = !(m_use_mac_phi_in_godunov);
    compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_new_const(),
                       get_density_new_const(), get_tracer_new_const(), get_tracer_new_const(),
                       include_pressure_gradient);
    compute_MAC_projected_velocities(get_velocity_new_const(), get_density_new_const(),
                                     AMREX_D_DECL(GetVecOfPtrs(u_mac), GetVecOfPtrs(v_mac),
                                     GetVecOfPtrs(w_mac)), GetVecOfPtrs(vel_forces), new_time);
    // *************************************************************************************
    // Compute the explicit "new" advective terms R_u^(n+1,*), R_r^(n+1,*) and R_t^(n+1,*)
    // *************************************************************************************
    compute_convective_term(get_conv_velocity_new(), get_conv_density_new(), get_conv_tracer_new(),
                            get_velocity_new_const(), get_density_new_const(), get_tracer_new_const(),
                            AMREX_D_DECL(GetVecOfPtrs(u_mac), GetVecOfPtrs(v_mac),
                            GetVecOfPtrs(w_mac)),
                            {}, {}, new_time);

    // *************************************************************************************
    // Compute viscosity / diffusive coefficients
    // *************************************************************************************
    compute_viscosity(GetVecOfPtrs(vel_eta), get_density_new(), get_velocity_new(), new_time, 1);

    // Here we create divtau of the (n+1,*) state that was computed in the predictor
    if ( (m_diff_type == DiffusionType::Explicit) || use_tensor_correction )
    {
        compute_divtau(get_divtau_new(), get_velocity_new_const(),
                       get_density_new_const(), GetVecOfConstPtrs(vel_eta));
    }

    // *************************************************************************************
    // Update density
    // *************************************************************************************
    update_density(StepType::Corrector);

    // *************************************************************************************
    // Update tracer
    // *************************************************************************************
    update_tracer(StepType::Corrector, tra_eta, tra_forces);

    // *************************************************************************************
    // Update velocity
    // *************************************************************************************
    update_velocity(StepType::Corrector, vel_eta, vel_forces);

    // **********************************************************************************************
    // Project velocity field, update pressure
    // **********************************************************************************************
    bool incremental_projection = false;
    ApplyProjection(get_density_nph_const(), new_time,m_dt,incremental_projection);

#ifdef AMREX_USE_EB
    // **********************************************************************************************
    //
    // Over-write velocity in cells with vfrac < 1e-4
    //
    // **********************************************************************************************
    incflo_correct_small_cells(get_velocity_new(),
                               AMREX_D_DECL(GetVecOfConstPtrs(u_mac), GetVecOfConstPtrs(v_mac),
                               GetVecOfConstPtrs(w_mac)));
#endif

#ifdef INCFLO_USE_PARTICLES
    // **************************************************************************************
    // Update the particle positions
    // **************************************************************************************
    evolveTracerParticles(AMREX_D_DECL(GetVecOfConstPtrs(u_mac), GetVecOfConstPtrs(v_mac),
                                       GetVecOfConstPtrs(w_mac)));
#endif
}
