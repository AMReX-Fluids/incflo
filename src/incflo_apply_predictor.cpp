#include <incflo.H>

using namespace amrex;
//
// Apply predictor:
//
//  1. Use u = vel_old to compute
//
//      if (!advect_momentum) then
//          conv_u  = - u grad u
//      else
//          conv_u  = - del dot (rho u u)
//      conv_r  = - div( u rho  )
//      if (m_iconserv_tracer) then
//          conv_t  = - div( u trac )
//      else
//          conv_t  = - u dot grad trac
//      eta_old     = visosity at m_cur_time
//      if (m_diff_type == DiffusionType::Explicit)
//         divtau _old = div( eta ( (grad u) + (grad u)^T ) ) / rho^n
//         rhs = u + dt * ( conv + divtau_old )
//      else
//         divtau_old  = 0.0
//         rhs = u + dt * conv
//
//      eta     = eta at new_time
//
//  2. Add explicit forcing term i.e. gravity + lagged pressure gradient
//
//      rhs += dt * ( g - grad(p + p0) / rho^nph )
//
//  3. A. If (m_diff_type == DiffusionType::Implicit)
//        solve implicit diffusion equation for u*
//
//     ( 1 - dt / rho^nph * div ( eta grad ) ) u* = u^n + dt * conv_u
//                                                  + dt * ( g - grad(p + p0) / rho^nph )
//
//     B. If (m_diff_type == DiffusionType::Crank-Nicolson)
//        solve semi-implicit diffusion equation for u*
//
//     ( 1 - (dt/2) / rho^nph * div ( eta_old grad ) ) u* = u^n +
//            dt * conv_u + (dt/2) / rho * div (eta_old grad) u^n
//          + dt * ( g - grad(p + p0) / rho^nph )
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
//     div( grad(phi) / rho^nph ) = div( u** )
//
//     Update pressure:
//
//     p = phi / dt
//
//     Update velocity, now divergence free
//
//     vel = u** - dt * grad p / rho^nph
//
// It is assumed that the ghost cells of the old data have been filled and
// the old and new data are the same in valid region.
//
void incflo::ApplyPredictor (bool incremental_projection)
{
    BL_PROFILE("incflo::ApplyPredictor");

    // We use the new time value for things computed on the "*" state
    Real new_time = m_cur_time + m_dt;

    // *************************************************************************************
    // Allocate space for the MAC velocities
    // *************************************************************************************
    Vector<MultiFab> u_mac(finest_level+1), v_mac(finest_level+1), w_mac(finest_level+1);
    int ngmac = nghost_mac();

    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(u_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(0)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev));,
                     v_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(1)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev));,
                     w_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(2)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev)););
        // do we still want to do this now that we always call a FillPatch (and all ghost cells get filled)?
        if (ngmac > 0) {
            AMREX_D_TERM(u_mac[lev].setBndry(0.0);,
                         v_mac[lev].setBndry(0.0);,
                         w_mac[lev].setBndry(0.0););
        }
    }

    // *************************************************************************************
    // Allocate space for half-time density
    // *************************************************************************************
    // Forcing terms
    Vector<MultiFab> vel_forces, tra_forces;

    Vector<MultiFab> vel_eta, tra_eta;

    // *************************************************************************************
    // Allocate space for the forcing terms
    // *************************************************************************************
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
    // Compute viscosity / diffusive coefficients
    // *************************************************************************************
    compute_viscosity(GetVecOfPtrs(vel_eta),
                      get_density_old(), get_velocity_old(),get_tracer_old(),
                      m_cur_time, 1);
    //when VOF method is used to advect the tracer, density and viscosity of each cell will
    // depend the VOF field value of the cell.
    if (m_vof_advect_tracer)
        update_vof_density (get_density_old(),get_tracer_old());

    // *************************************************************************************
    // Compute explicit viscous term
    // Note that for !advect_momentum, this actually computes divtau / rho
    // *************************************************************************************
    if (need_divtau() || use_tensor_correction )
    {
        compute_divtau(get_divtau_old(),get_velocity_old_const(),
                       get_density_old_const(),GetVecOfConstPtrs(vel_eta));
    }

    // *************************************************************************************
    // Compute explicit diffusive term -- note this is used inside compute_convective_term
    // *************************************************************************************
    if (m_advect_tracer)
    {
        compute_tracer_diff_coeff(GetVecOfPtrs(tra_eta),1);
        if (need_divtau()) {
            compute_laps(get_laps_old(), get_tracer_old_const(), GetVecOfConstPtrs(tra_eta));
        }
    }

    // **********************************************************************************************
    // Compute the forcing terms
    // *************************************************************************************
    bool include_pressure_gradient = !(m_use_mac_phi_in_godunov);
    compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_old_const(),
                       get_density_old_const(), get_tracer_old_const(), get_tracer_old_const(),
                       include_pressure_gradient);

    // **********************************************************************************************
    // Compute the MAC-projected velocities at all levels
    // *************************************************************************************
    compute_MAC_projected_velocities(get_velocity_old_const(), get_density_old_const(),
                                     AMREX_D_DECL(GetVecOfPtrs(u_mac), GetVecOfPtrs(v_mac),
                                     GetVecOfPtrs(w_mac)), GetVecOfPtrs(vel_forces), m_cur_time);

    // *************************************************************************************
    // if (advection_type == "Godunov")
    //      Compute the explicit advective terms R_u^(n+1/2), R_s^(n+1/2) and R_t^(n+1/2)
    // if (advection_type == "MOL"                )
    //      Compute the explicit advective terms R_u^n      , R_s^n       and R_t^n
    // Note that if advection_type != "MOL" then we call compute_tra_forces inside this routine
    // *************************************************************************************
    compute_convective_term(get_conv_velocity_old(), get_conv_density_old(), get_conv_tracer_old(),
                            get_velocity_old_const(), get_density_old_const(), get_tracer_old_const(),
                            AMREX_D_DECL(GetVecOfPtrs(u_mac), GetVecOfPtrs(v_mac),
                            GetVecOfPtrs(w_mac)),
                            GetVecOfPtrs(vel_forces), GetVecOfPtrs(tra_forces),
                            m_cur_time);

    // *************************************************************************************
    // Update density
    // *************************************************************************************
    update_density(StepType::Predictor);

    // **********************************************************************************************
    // Update tracer
    // **********************************************************************************************
    update_tracer(StepType::Predictor, tra_eta, tra_forces);

    // **********************************************************************************************
    // Update velocity
    // **********************************************************************************************
    update_velocity(StepType::Predictor, vel_eta, vel_forces);

    // **********************************************************************************************
    // Project velocity field, update pressure
    // **********************************************************************************************
    ApplyProjection(get_density_nph_const(),
                    AMREX_D_DECL(GetVecOfPtrs(u_mac), GetVecOfPtrs(v_mac),
                    GetVecOfPtrs(w_mac)),new_time,m_dt,incremental_projection);

#ifdef INCFLO_USE_PARTICLES
    // **************************************************************************************
    // Update the particle positions
    // **************************************************************************************
    if (m_advection_type != "MOL") {
        evolveTracerParticles(AMREX_D_DECL(GetVecOfConstPtrs(u_mac), GetVecOfConstPtrs(v_mac),
                                           GetVecOfConstPtrs(w_mac)));
    }
#endif

#ifdef AMREX_USE_EB
    // **********************************************************************************************
    // Over-write velocity in cells with vfrac < 1e-4
    // **********************************************************************************************
    if (m_advection_type == "MOL")
        incflo_correct_small_cells(get_velocity_new(),
                                   AMREX_D_DECL(GetVecOfConstPtrs(u_mac), GetVecOfConstPtrs(v_mac),
                                   GetVecOfConstPtrs(w_mac)));
#endif

// use vof to advect tracer
    if (!incremental_projection)
      tracer_vof_advection(get_tracer_new (), AMREX_D_DECL(GetVecOfConstPtrs(u_mac), GetVecOfConstPtrs(v_mac),
                           GetVecOfConstPtrs(w_mac)));

}
