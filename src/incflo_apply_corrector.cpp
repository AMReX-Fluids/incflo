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

    constexpr Real m_half = Real(0.5);

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

    // *************************************************************************************
    // Allocate space for half-time density
    // *************************************************************************************
    Vector<MultiFab> density_nph;
    for (int lev = 0; lev <= finest_level; ++lev) {
        density_nph.emplace_back(grids[lev], dmap[lev], 1, 0, MFInfo(), Factory(lev));
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

    // **********************************************************************************************
    // Compute the MAC-projected velocities at all levels
    // *************************************************************************************
    bool include_pressure_gradient = !(m_use_mac_phi_in_godunov);
    compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_new_const(),
                       get_density_new_const(), get_tracer_new_const(), get_tracer_new_const(),
                       include_pressure_gradient);
    compute_MAC_projected_velocities(get_velocity_new_const(), get_density_new_const(),
                                     AMREX_D_DECL(GetVecOfPtrs(u_mac), GetVecOfPtrs(v_mac),
                                     GetVecOfPtrs(w_mac)), GetVecOfPtrs(vel_forces), new_time);
    // **********************************************************************************************
    // Compute the explicit "new" advective terms R_u^(n+1,*), R_r^(n+1,*) and R_t^(n+1,*)
    // Note that "get_conv_tracer_new" returns div(rho u tracer)
    // *************************************************************************************
    compute_convective_term(get_conv_velocity_new(), get_conv_density_new(), get_conv_tracer_new(),
                            get_velocity_new_const(), get_density_new_const(), get_tracer_new_const(),
                            AMREX_D_DECL(GetVecOfPtrs(u_mac), GetVecOfPtrs(v_mac),
                            GetVecOfPtrs(w_mac)),
                            {}, {}, new_time);

    // *************************************************************************************
    // Compute viscosity / diffusive coefficients
    // *************************************************************************************
    compute_viscosity(GetVecOfPtrs(vel_eta),
                      get_density_new(), get_velocity_new(),
                      new_time, 1);
    compute_tracer_diff_coeff(GetVecOfPtrs(tra_eta),1);

    // Here we create divtau of the (n+1,*) state that was computed in the predictor;
    //      we use this laps only if DiffusionType::Explicit
    if ( (m_diff_type == DiffusionType::Explicit) || use_tensor_correction )
    {
        compute_divtau(get_divtau_new(), get_velocity_new_const(),
                       get_density_new_const(), GetVecOfConstPtrs(vel_eta));
    }

    if (m_advect_tracer && m_diff_type == DiffusionType::Explicit) {
        compute_laps(get_laps_new(), get_tracer_new_const(), GetVecOfConstPtrs(tra_eta));
    }

    // *************************************************************************************
    // Define local variables for lambda to capture.
    // *************************************************************************************
    Real l_dt = m_dt;
    int l_ntrac = (m_advect_tracer) ? m_ntrac : 0;

    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (m_constant_density)
    {
        for (int lev = 0; lev <= finest_level; lev++)
            MultiFab::Copy(density_nph[lev], m_leveldata[lev]->density_o, 0, 0, 1, 0);
    } else {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real const> const& rho_o  = ld.density_o.const_array(mfi);
                Array4<Real> const& rho_n        = ld.density.array(mfi);
                Array4<Real> const& rho_nph      = density_nph[lev].array(mfi);
                Array4<Real const> const& drdt_o = ld.conv_density_o.const_array(mfi);
                Array4<Real const> const& drdt   = ld.conv_density.const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const Real rho_old = rho_o(i,j,k);

                    Real rho_new = rho_old + l_dt * m_half*(drdt(i,j,k)+drdt_o(i,j,k));
                    rho_nph(i,j,k) = m_half * (rho_old + rho_new);

                    rho_n  (i,j,k) = rho_new;
                });
            } // mfi
        } // lev

        // Average down solution
        for (int lev = finest_level-1; lev >= 0; --lev) {
#ifdef AMREX_USE_EB
            amrex::EB_average_down(m_leveldata[lev+1]->density, m_leveldata[lev]->density,
                                   0, 1, refRatio(lev));
#else
            amrex::average_down(m_leveldata[lev+1]->density, m_leveldata[lev]->density,
                                0, 1, refRatio(lev));
#endif
        }
    } // not constant density

    // *************************************************************************************
    // Compute the tracer forcing terms (forcing for (rho s), not for s)
    // *************************************************************************************
    if (m_advect_tracer)
        compute_tra_forces(GetVecOfPtrs(tra_forces),  GetVecOfConstPtrs(density_nph));

    // *************************************************************************************
    // Update the tracer next (note that dtdt already has rho in it)
    // (rho trac)^new = (rho trac)^old + dt * (
    //                   div(rho trac u) + div (mu grad trac) + rho * f_t
    // *************************************************************************************
    if (m_advect_tracer)
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real const> const& tra_o   = ld.tracer_o.const_array(mfi);
                Array4<Real const> const& rho_o   = ld.density_o.const_array(mfi);
                Array4<Real      > const& tra     = ld.tracer.array(mfi);
                Array4<Real const> const& rho     = ld.density.const_array(mfi);
                Array4<Real const> const& dtdt_o  = ld.conv_tracer_o.const_array(mfi);
                Array4<Real const> const& dtdt    = ld.conv_tracer.const_array(mfi);
                Array4<Real const> const& tra_f   = (l_ntrac > 0) ? tra_forces[lev].const_array(mfi)
                                                                : Array4<Real const>{};
                auto iconserv = get_tracer_iconserv_device_ptr();

                if (m_diff_type == DiffusionType::Explicit)
                {
                    Array4<Real const> const& laps_o = (l_ntrac > 0) ? ld.laps_o.const_array(mfi)
                                                                     : Array4<Real const>{};
                    Array4<Real const> const& laps   = (l_ntrac > 0) ? ld.laps.const_array(mfi)
                                                                     : Array4<Real const>{};
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < l_ntrac; ++n)
                        {
                            if ( iconserv[n] ) {
                                Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt * (
                                     m_half*(  dtdt(i,j,k,n) + dtdt_o(i,j,k,n))
                                    +m_half*(laps_o(i,j,k,n) +   laps(i,j,k,n))
                                    +    tra_f(i,j,k,n) );

                                tra(i,j,k,n) = tra_new / rho(i,j,k);
                            } else {
                                tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt * (
                                     m_half*(  dtdt(i,j,k,n) + dtdt_o(i,j,k,n))
                                    +m_half*(laps_o(i,j,k,n) +   laps(i,j,k,n))
                                    +    tra_f(i,j,k,n) );
                            }
                        }
                    });
                }
                else if (m_diff_type == DiffusionType::Crank_Nicolson)
                {
                    Array4<Real const> const& laps_o = (l_ntrac > 0) ? ld.laps_o.const_array(mfi)
                                                                     : Array4<Real const>{};
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < l_ntrac; ++n)
                        {
                            if ( iconserv[n] ) {
                                Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt * (
                                     m_half*(  dtdt(i,j,k,n) + dtdt_o(i,j,k,n))
                                    +m_half*(laps_o(i,j,k,n)                  )
                                    +    tra_f(i,j,k,n) );

                                tra(i,j,k,n) = tra_new / rho(i,j,k);
                            } else {
                                tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt * (
                                     m_half*(  dtdt(i,j,k,n) + dtdt_o(i,j,k,n))
                                    +m_half*(laps_o(i,j,k,n)                  )
                                    +    tra_f(i,j,k,n) );
                            }
                    });
                }
                else if (m_diff_type == DiffusionType::Implicit)
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < l_ntrac; ++n)
                        {
                            if ( iconserv[n] ) {
                                Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt * (
                                     m_half*(  dtdt(i,j,k,n)+dtdt_o(i,j,k,n))
                                    +      tra_f(i,j,k,n) );

                                tra(i,j,k,n) = tra_new / rho(i,j,k);
                            } else {
                                tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt * (
                                     m_half*(  dtdt(i,j,k,n)+dtdt_o(i,j,k,n))
                                    +      tra_f(i,j,k,n) );
                            }
                    });
                }
            } // mfi
        } // lev
    } // if (m_advect_tracer)

    // *************************************************************************************
    // Solve diffusion equation for tracer
    // *************************************************************************************
    if ( m_advect_tracer )
    {
        if (m_diff_type == DiffusionType::Crank_Nicolson || m_diff_type == DiffusionType::Implicit)
        {
            const int ng_diffusion = 1;
            for (int lev = 0; lev <= finest_level; ++lev)
                fillphysbc_tracer(lev, new_time, m_leveldata[lev]->tracer, ng_diffusion);

            Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_dt : m_half*m_dt;
            diffuse_scalar(get_tracer_new(), get_density_new(), GetVecOfConstPtrs(tra_eta), dt_diff);
        }
        else
        {
            // Need to average down tracer since the diffusion solver didn't do it for us.
            for (int lev = finest_level-1; lev >= 0; --lev) {
#ifdef AMREX_USE_EB
                amrex::EB_average_down(m_leveldata[lev+1]->tracer, m_leveldata[lev]->tracer,
                                       0, m_ntrac, refRatio(lev));
#else
                amrex::average_down(m_leveldata[lev+1]->tracer, m_leveldata[lev]->tracer,
                                    0, m_ntrac, refRatio(lev));
#endif
            }
        }
    } // if (m_advect_tracer)


    // *************************************************************************************
    // Define the forcing terms to use in the final update (using half-time density)
    // *************************************************************************************
    compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_new_const(),
                       GetVecOfConstPtrs(density_nph),
                       get_tracer_old_const(), get_tracer_new_const());

    // *************************************************************************************
    // Update velocity
    // *************************************************************************************
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& vel = ld.velocity.array(mfi);
            Array4<Real const> const& vel_o = ld.velocity_o.const_array(mfi);
            Array4<Real const> const& dvdt = ld.conv_velocity.const_array(mfi);
            Array4<Real const> const& dvdt_o = ld.conv_velocity_o.const_array(mfi);
            Array4<Real const> const& vel_f = vel_forces[lev].const_array(mfi);

            Array4<Real const> const& rho_old  = ld.density_o.const_array(mfi);
            Array4<Real const> const& rho_new  = ld.density.const_array(mfi);
            Array4<Real const> const& rho_nph  = density_nph[lev].array(mfi);

            if (m_diff_type == DiffusionType::Explicit)
            {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                Array4<Real const> const& divtau   = ld.divtau.const_array(mfi);

                if (m_advect_momentum) {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                    m_half*(  dvdt_o(i,j,k,0)+  dvdt(i,j,k,0))
                                                  + m_half*(divtau_o(i,j,k,0)+divtau(i,j,k,0))
                                                  + rho_nph(i,j,k) * vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) =  rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                     m_half*(  dvdt_o(i,j,k,1)+  dvdt(i,j,k,1))
                                                   + m_half*(divtau_o(i,j,k,1)+divtau(i,j,k,1))
                                                   + rho_nph(i,j,k) * vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) =  rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                     m_half*(  dvdt_o(i,j,k,2)+  dvdt(i,j,k,2))
                                                   + m_half*(divtau_o(i,j,k,2)+divtau(i,j,k,2))
                                                   + rho_nph(i,j,k) * vel_f(i,j,k,2) ););

                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                    m_half*(  dvdt_o(i,j,k,0)+  dvdt(i,j,k,0))
                                                  + m_half*(divtau_o(i,j,k,0)+divtau(i,j,k,0))
                                                  + vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                    m_half*(  dvdt_o(i,j,k,1)+  dvdt(i,j,k,1))
                                                  + m_half*(divtau_o(i,j,k,1)+divtau(i,j,k,1))
                                                  + vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                    m_half*(  dvdt_o(i,j,k,2)+  dvdt(i,j,k,2))
                                                  + m_half*(divtau_o(i,j,k,2)+divtau(i,j,k,2))
                                                  + vel_f(i,j,k,2) ););
                    });
                }
            }
            else if (m_diff_type == DiffusionType::Crank_Nicolson)
            {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);

                if (m_advect_momentum) {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                    m_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0))
                                                   +m_half*(divtau_o(i,j,k,0) ) + rho_nph(i,j,k)*vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) = rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                    m_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1))
                                                   +m_half*(divtau_o(i,j,k,1) ) + rho_nph(i,j,k)*vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) = rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                    m_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2))
                                                   +m_half*(divtau_o(i,j,k,2) ) + rho_nph(i,j,k)*vel_f(i,j,k,2) ););

                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                    m_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0))
                                                  + m_half* divtau_o(i,j,k,0) + vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                    m_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1))
                                                  + m_half* divtau_o(i,j,k,1) + vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                    m_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2))
                                                  + m_half* divtau_o(i,j,k,2) + vel_f(i,j,k,2) ););
                    });
                }
            }
            else if (m_diff_type == DiffusionType::Implicit)
            {
                if (use_tensor_correction)
                {
                    Array4<Real const> const& divtau   = ld.divtau.const_array(mfi);
                    if (m_advect_momentum)
                    {
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                        m_half*(dvdt_o(i,j,k,0) +   dvdt(i,j,k,0))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,0) + divtau(i,j,k,0));,
//
                                         vel(i,j,k,1) = rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                        m_half*(dvdt_o(i,j,k,1) +   dvdt(i,j,k,1))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,1) + divtau(i,j,k,1));,
//
                                         vel(i,j,k,2) = rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                        m_half*(dvdt_o(i,j,k,2) +   dvdt(i,j,k,2))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,2) + divtau(i,j,k,2)););

                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                        m_half*(  dvdt_o(i,j,k,0) +   dvdt(i,j,k,0))
                                                      + vel_f(i,j,k,0) + divtau(i,j,k,0));,
//
                                         vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                        m_half*(dvdt_o(i,j,k,1) +   dvdt(i,j,k,1))
                                                      + vel_f(i,j,k,1) + divtau(i,j,k,1) );,
//
                                         vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                        m_half*(  dvdt_o(i,j,k,2) +   dvdt(i,j,k,2))
                                                      + vel_f(i,j,k,2) + divtau(i,j,k,2) ););
                        });
                    }
                } else {
                    if (m_advect_momentum)
                    {
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                        m_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,0) );,
                                         vel(i,j,k,1) = rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                        m_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,1) );,
                                         vel(i,j,k,2) = rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                        m_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,2) ););
                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                        m_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0)) + vel_f(i,j,k,0) );,
                                         vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                        m_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1)) + vel_f(i,j,k,1) );,
                                         vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                        m_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2)) + vel_f(i,j,k,2) ););
                        });
                    }
                }
            }
        }
    }

    // **********************************************************************************************
    //
    // Solve diffusion equation for u* at t^{n+1} but using eta at predicted new time
    //
    // **********************************************************************************************

    if (m_diff_type == DiffusionType::Crank_Nicolson || m_diff_type == DiffusionType::Implicit)
    {
        const int ng_diffusion = 1;
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            fillphysbc_velocity(lev, new_time, m_leveldata[lev]->velocity, ng_diffusion);
            fillphysbc_density (lev, new_time, m_leveldata[lev]->density , ng_diffusion);
        }

        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_dt : m_half*m_dt;
        diffuse_velocity(get_velocity_new(), get_density_new(), GetVecOfConstPtrs(vel_eta), dt_diff);
    }

    // **********************************************************************************************
    //
    // Project velocity field, update pressure
    bool incremental_projection = false;
    ApplyProjection(GetVecOfConstPtrs(density_nph), new_time,m_dt,incremental_projection);

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
}
