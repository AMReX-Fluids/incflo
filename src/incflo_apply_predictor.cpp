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
//      conv_t  = - div( u trac )
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
// It is assumed that the ghost cels of the old data have been filled and
// the old and new data are the same in valid region.
//
void incflo::ApplyPredictor (bool incremental_projection)
{
    BL_PROFILE("incflo::ApplyPredictor");

    constexpr Real m_half = Real(0.5);

    // We use the new time value for things computed on the "*" state
    Real new_time = m_cur_time + m_dt;

#ifdef INCFLO_USE_MOVING_EB
    // FIXME - would these be good from the corrector step? can we rely on doing a corrector?
    //   what about initial iterations
    // *************************************************************************************
    // Reset the solvers to work with the new EB
    // *************************************************************************************
    m_diffusion_tensor_op.reset();
    m_diffusion_scalar_op.reset();

    macproj.reset(new Hydro::MacProjector(Geom(0,finest_level),
                      MLMG::Location::FaceCentroid,  // Location of mac_vec
                      MLMG::Location::FaceCentroid,  // Location of beta
                      MLMG::Location::CellCenter  ) ); // Location of solution variable phi
#endif

    if (m_verbose > 2)
    {
        amrex::Print() << "Before predictor step:" << std::endl;
        //PrintMaxValues(new_time);
    }

    // *************************************************************************************
    // Allocate space for the MAC velocities
    // *************************************************************************************
    Vector<MultiFab> u_mac(finest_level+1), v_mac(finest_level+1), w_mac(finest_level+1);
    int ngmac = nghost_mac();

    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(u_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(0)), dmap[lev],
                          1, ngmac, MFInfo(), OldFactory(lev));,
                     v_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(1)), dmap[lev],
                          1, ngmac, MFInfo(), OldFactory(lev));,
                     w_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(2)), dmap[lev],
                          1, ngmac, MFInfo(), OldFactory(lev)););
        // do we still want to do this now that we always call a FillPatch (and all ghost cells get filled)?
        if (ngmac > 0) {
            AMREX_D_TERM(u_mac[lev].setBndry(0.0);,
                         v_mac[lev].setBndry(0.0);,
                         w_mac[lev].setBndry(0.0););
        }
    }
    // *************************************************************************************
    // Allocate space for half-time density using "old" EB
    // *************************************************************************************
    Vector<MultiFab> density_nph_oldeb;
    for (int lev = 0; lev <= finest_level; ++lev)
        density_nph_oldeb.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(), OldFactory(lev));

    // *************************************************************************************
    // Forcing terms
    // *************************************************************************************
    Vector<MultiFab> vel_forces, tra_forces;

    Vector<MultiFab> vel_eta, tra_eta;

    // *************************************************************************************
    // Allocate space for the forcing terms
    // *************************************************************************************
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel_forces.emplace_back(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost_force(),
                                MFInfo(), OldFactory(lev));

        if (m_advect_tracer) {
            tra_forces.emplace_back(grids[lev], dmap[lev], m_ntrac, nghost_force(),
                                    MFInfo(), OldFactory(lev));
        }
        vel_eta.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(), OldFactory(lev));
        if (m_advect_tracer) {
            tra_eta.emplace_back(grids[lev], dmap[lev], m_ntrac, 1, MFInfo(), OldFactory(lev));
        }
    }

    // *************************************************************************************
    // We now define the forcing terms to use in the Godunov prediction inside the predictor
    // *************************************************************************************

    // *************************************************************************************
    // Compute viscosity / diffusive coefficients
    // *************************************************************************************
    compute_viscosity(GetVecOfPtrs(vel_eta),
                      get_density_old(), get_velocity_old(),
                      m_cur_time, 1);
    compute_tracer_diff_coeff(GetVecOfPtrs(tra_eta),1);

    // *************************************************************************************
    // Compute explicit viscous term
    // *************************************************************************************
    if (need_divtau() || use_tensor_correction)
    {
        compute_divtau(get_divtau_old(),get_velocity_old_const(),
                       get_density_old_const(),GetVecOfConstPtrs(vel_eta));
    }

    // *************************************************************************************
    // Compute explicit diffusive terms
    // *************************************************************************************
    if (m_advect_tracer && need_divtau())
    {
        compute_laps(get_laps_old(), get_tracer_old_const(), get_density_old_const(),
                     GetVecOfConstPtrs(tra_eta));
    }

    // *************************************************************************************
    // Compute the MAC-projected velocities at all levels
    // *************************************************************************************
    bool include_pressure_gradient = !(m_use_mac_phi_in_godunov);
    compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_old_const(),
                       get_density_old_const(), get_tracer_old_const(), get_tracer_new_const(),
                       include_pressure_gradient);

    //VisMF::Write(m_leveldata[0]->density_o,"do3");
    //amrex::Print() << "density: " << m_leveldata[0]->density_o << std::endl;

    compute_MAC_projected_velocities(get_velocity_old_const(), get_density_old_const(),
                                     AMREX_D_DECL(GetVecOfPtrs(u_mac), GetVecOfPtrs(v_mac),
                                     GetVecOfPtrs(w_mac)), GetVecOfPtrs(vel_forces), m_cur_time);

    //FIXME this density is not the same as the print out from first thing inside compute_convective_term
    //VisMF::Write(m_leveldata[0]->density_o,"do4");
    //amrex::Print() << "desnsity_o:\n" << m_leveldata[0]->density_o << std::endl;

    // *************************************************************************************
    // if (advection_type == "Godunov")
    //      Compute the explicit advective terms R_u^(n+1/2), R_s^(n+1/2) and R_t^(n+1/2)
    // if (advection_type == "MOL"                )
    //      Compute the explicit advective terms R_u^n      , R_s^n       and R_t^n
    // Note that "get_conv_tracer_old" returns div(rho u tracer)
    // *************************************************************************************

    compute_convective_term(get_conv_velocity_old(), get_conv_density_old(), get_conv_tracer_old(),
                            get_velocity_old_const(), get_density_old_const(), get_tracer_old_const(),
                            AMREX_D_DECL(GetVecOfPtrs(u_mac), GetVecOfPtrs(v_mac),
                            GetVecOfPtrs(w_mac)),
                            GetVecOfPtrs(vel_forces), GetVecOfPtrs(tra_forces),
                            m_cur_time);

    // *************************************************************************************
    // Define local variables for lambda to capture.
    // *************************************************************************************
    Real l_dt = m_dt;
    bool l_constant_density = m_constant_density;

    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (l_constant_density)
    {
        for (int lev = 0; lev <= finest_level; lev++)
            MultiFab::Copy(density_nph_oldeb[lev], m_leveldata[lev]->density_o, 0, 0, 1, 1);
    } else {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];

#ifdef AMREX_USE_MOVING_EB
            //
	    // For moving EB, redistribute and returns full state at new time
            //
	    redistribute_term(ld.density, ld.conv_density_o, ld.density_o,
			      get_density_bcrec_device_ptr(), lev,
			      get_velocity_eb()[lev]);
#else
	    //
	    // Add advective update that's already been redistributed to get rho^(n+1)
	    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real  const> const& rho_o  = ld.density_o.const_array(mfi);
                Array4<Real> const& rho_new       = ld.density.array(mfi);
                Array4<Real const> const& drdt    = ld.conv_density_o.const_array(mfi);

                auto const& vfrac_old = OldEBFactory(lev).getVolFrac().const_array(mfi);
                auto const& vfrac_new =    EBFactory(lev).getVolFrac().const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rho_new(i,j,k) =  rho_o(i,j,k) + l_dt * drdt(i,j,k);
		});
            } // mfi
#endif

            // Fill ghost cells of the new density field so that we can define density_nph
            //      on the valid region grown by 1
            int ng = 1;
            fillpatch_density(lev, m_t_new[lev], ld.density, ng);


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                 Box const& bxg1 = mfi.growntilebox(1);
                 Array4<Real  const> const& rho_old  = ld.density_o.const_array(mfi);
                 Array4<Real  const> const& rho_new  = ld.density.const_array(mfi);
                 Array4<Real>        const& rho_nph  = density_nph_oldeb[lev].array(mfi);

                 auto const& vfrac_new = OldEBFactory(lev).getVolFrac().const_array(mfi);
                 auto const& vfrac_old =    EBFactory(lev).getVolFrac().const_array(mfi);

                 amrex::ParallelFor(bxg1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                 {
                     // FIXME - probably want to think more about what to do with rho_n+1/2 here
                     // if (vfrac_new(i,j,k) > 0. && vfrac_old(i,j,k) == 0.)
                     // {
                     //          rho_nph(i,j,k) = rho_new(i,j,k);
                     // } else {
                         rho_nph(i,j,k) = m_half * (rho_old(i,j,k) + rho_new(i,j,k));
                     // }
                     // if ((vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1.))
                     //     amrex::Print() << "rho" << IntVect(i,j) << ": " << rho_new(i,j,0) << std::endl;
                 });

#ifdef AMREX_USE_MOVING_EB
            //
	    // Create S^m(rho) and update-hat
	    //
                Box const& bx = mfi.tilebox();
		// Store update-hat here...
		Array4<Real      > const& drdt_o   = ld.conv_density_o.array(mfi);

		Box const& gbx = amrex::grow(bx,ld.density.nGrow());
		FArrayBox update_fab(gbx,1,The_Async_Arena());
		Array4<Real> const& update = update_fab.array();
		amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		{
		    update(i,j,k) = 0.;
		});

		// redistribute - own lambda...
		redistribute_term(mfi, drdt_o, update, rho_old,
				  get_density_bcrec_device_ptr(), lev,
				  get_velocity_eb()[lev]->const_array(mfi));

		// Make update-hat
		// NOt saving S^m(rho) for now, don't think we need it
		// FIXME - I think this is okay as redistribute only uses tra on bx
		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		{
		    // Update-hat = S^m(rho^n) - S^m(rho*trac^n - dt A^n)
		    //            = S^m(rho^n) - rho^(p,n+1)
		    drdt_o(i,j,k) -= rho_new(i,j,k);
		    drdt_o(i,j,k) /= l_dt;
		});
#endif
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
    // Compute (or if Godunov, re-compute) the tracer forcing terms (forcing for (rho s), not for s)
    // *************************************************************************************
    if (m_advect_tracer)
       compute_tra_forces(GetVecOfPtrs(tra_forces), GetVecOfConstPtrs(density_nph_oldeb));

    // *************************************************************************************
    // Update the tracer next
    // *************************************************************************************
    int l_ntrac = (m_advect_tracer) ? m_ntrac : 0;

    if (m_advect_tracer)
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];

#ifdef AMREX_USE_MOVING_EB
	    //
	    // For moving EB, assemble state for redistribution
	    //
	    MultiFab update(grids[lev], dmap[lev], m_ntrac, nghost_state(), MFInfo(), Factory(lev));

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& update_t     = update.array(mfi);
                Array4<Real const> const& dtdt_o = ld.conv_tracer_o.const_array(mfi);
// FIXME - why do we have this test when init.cpp requires at least 1 tracer...
		Array4<Real const> const& laps_o = (l_ntrac > 0) ? ld.laps_o.const_array(mfi)
                                                                 : Array4<Real const>{};
                Array4<Real const> const& tra_f   = (l_ntrac > 0) ? tra_forces[lev].const_array(mfi)
                                                                : Array4<Real const>{};

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
		    // Build update
		    for ( int n = 0; n < l_ntrac; n++)
		    {
			update_t(i,j,k,n) = dtdt_o(i,j,k,n) + tra_f(i,j,k,n) + laps_o(i,j,k,n);
		    }
		});
	    }

	    // Need to ensure that boundaries are consistent
	    update.FillBoundary(geom[lev].periodicity());
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
		Array4<Real> const& tra          = ld.tracer.array(mfi);
                Array4<Real const> const& rho    = ld.density.const_array(mfi);
                Array4<Real const> const& rho_o  = ld.density_o.const_array(mfi);
		Array4<Real const> const& dtdt_o  = ld.conv_tracer_o.const_array(mfi);
		//FIXME? - init will abort if ntrac < 1...
                Array4<Real const> const& tra_f  = (l_ntrac > 0) ? tra_forces[lev].const_array(mfi)
		                                                 : Array4<Real const>{};

#ifdef AMREX_USE_MOVING_EB
		Array4<Real      > const& tra_o  = ld.tracer_o.array(mfi);
		Array4<Real const> const& laps_o = (l_ntrac > 0) ? ld.laps_o.const_array(mfi)
                                                                 : Array4<Real const>{};
                Array4<Real> const& update_t     = update.array(mfi);
		
                // Redistribute
		// (rho trac)^new = (rho trac)^old + dt * (
		//                   div(rho trac u) + div (mu grad trac) + rho * f_t

		// FIXME - I need to grow this box -- maybe safer to have 3 mfiters???
		// incflo_compute_advective... does this same thing...
		// except update now needs to be temporary, can't overwrite conv
		Box const& gbx = amrex::grow(bx,ld.tracer.nGrow());
		amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		{
		    for (int n = 0; n < l_ntrac; ++n)
		    {
			tra_o(i,j,k,n) = rho_o(i,j,k)*tra_o(i,j,k,n);
		    }
		});

		// redistribute - own lambda...
		redistribute_term(mfi, tra, update_t, tra_o, get_tracer_bcrec_device_ptr(), lev,
				  get_velocity_eb()[lev]->const_array(mfi));

                // Don't divide trac/rho yet..
		
#else

                Array4<Real const> const& tra_o   = ld.tracer_o.const_array(mfi);

                if (m_diff_type == DiffusionType::Explicit)
                {
                    Array4<Real const> const& laps_o = (l_ntrac > 0) ? ld.laps_o.const_array(mfi)
                                                                     : Array4<Real const>{};
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        // (rho trac)^new = (rho trac)^old + dt * (
                        //                   div(rho trac u) + div (mu grad trac) + rho * f_t
                        for (int n = 0; n < l_ntrac; ++n)
                        {
                            Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) + laps_o(i,j,k,n) );

                            tra_new /= rho(i,j,k);
                            tra(i,j,k,n) = tra_new;
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
                            Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) + m_half * laps_o(i,j,k,n) );

                            tra_new /= rho(i,j,k);
                            tra(i,j,k,n) = tra_new;
                        }
                    });
                }
                else if (m_diff_type == DiffusionType::Implicit)
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < l_ntrac; ++n)
                        {
                            Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) );

                            tra_new /= rho(i,j,k);
                            tra(i,j,k,n) = tra_new;
                        }
                    });
                }
#endif
            } // mfi

#ifdef AMREX_USE_MOVING_EB
            //
	    // Create S^m(rho*trac) and update-hat
	    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();

		Array4<Real      > const& rhotra_o = ld.tracer_o.array(mfi);
		Array4<Real      > const& rhotra   = ld.tracer.array(mfi);
		// Store update-hat here...
		Array4<Real      > const& dtdt_o   = ld.conv_tracer_o.array(mfi);

		// Why not n in parallel for?
		// FIXME - I need to grow this box -- maybe safer to have 3 mfiters???
		// incflo_compute_advective... does this same thing...
                // Could probably merge this MFIter with above as long as we make
		// temproary MFs for redistributed update-hat
		// but then we'd still need a copy, so probably doesn't matter
		Box const& gbx = amrex::grow(bx,ld.tracer.nGrow());
		FArrayBox update_fab(gbx,l_ntrac,The_Async_Arena());
		Array4<Real> const& zero_arr = update_fab.array();
		amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		{
		    for (int n = 0; n < l_ntrac; ++n)
		    {
			zero_arr(i,j,k,n) = 0.;
		    }
		});

		// redistribute - own lambda...
		redistribute_term(mfi, dtdt_o, zero_arr, rhotra_o,
				  get_tracer_bcrec_device_ptr(), lev,
				  get_velocity_eb()[lev]->const_array(mfi));

		// Make update-hat
		// NOt saving S^m(rho*tra) for now, don't think we need it
		// FIXME - I think this is okay as redistribute only uses tra on bx
		Array4<Real const> const& rho      = ld.density.const_array(mfi);
		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		{
		    for (int n = 0; n < l_ntrac; ++n)
		    {
			// Update-hat = S^m(rho*trac^n) - S^m(rho*trac^n - dt(A^n - D^n - f^n))
			//            = S^m(rho*trac^n) - rhotrac^(p,n+1)
			dtdt_o(i,j,k,n) -= rhotra(i,j,k,n);
			dtdt_o(i,j,k,n) /= l_dt;
			
			// put back trac = rhotrac/rho_new - non-ideal for MEB, but convection
			// expects just trac...
			rhotra(i,j,k,n) /= rho(i,j,k); 
		    }
		});
	    } // mfi
#endif
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
    // Define (or if advection_type != "MOL", re-define) the forcing terms, without the viscous terms
    //    and using the half-time density
    // *************************************************************************************
    compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_old_const(),
                       GetVecOfConstPtrs(density_nph_oldeb),
                       get_tracer_old_const(), get_tracer_new_const());


    // *************************************************************************************
    // Update the velocity
    // *************************************************************************************
    for (int lev = 0; lev <= finest_level; lev++)
    {
        auto& ld = *m_leveldata[lev];

#ifdef AMREX_USE_MOVING_EB
	//
	// For moving EB, assemble update for redistribution
	//
	// FIXME -- think we only need 3 here, not ng_state...
	MultiFab update(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost_state(),
			MFInfo(), Factory(lev));

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	    Box const& bx = mfi.tilebox();
	    Array4<Real> const& update_v       = update.array(mfi);
	    Array4<Real const> const& dvdt_o   = ld.conv_velocity_o.const_array(mfi);
	    Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
            Array4<Real const> const& vel_f    = vel_forces[lev].const_array(mfi);
	    Array4<Real const> const& rho_o    = ld.density_o.const_array(mfi);
	    
	    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		// Build update
		for ( int n = 0; n < AMREX_SPACEDIM; n++)
		{
		    update_v(i,j,k,n) = dvdt_o(i,j,k,n) + divtau_o(i,j,k,n)
			                + rho_o(i,j,k)*vel_f(i,j,k,n);
		}
	    });
	}

	// Need to ensure that boundaries are consistent
	update.FillBoundary(geom[lev].periodicity());
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real      > const& vel = ld.velocity.array(mfi);
            Array4<Real const> const& dvdt = ld.conv_velocity_o.const_array(mfi);
            Array4<Real const> const& vel_f = vel_forces[lev].const_array(mfi);
            Array4<Real const> const& rho_old  = ld.density_o.const_array(mfi);
            Array4<Real const> const& rho_new  = ld.density.const_array(mfi);
            Array4<Real const> const& rho_nph  = density_nph_oldeb[lev].array(mfi);

            // FIXME for debugging
            auto const& vfrac_old = OldEBFactory(lev).getVolFrac().const_array(mfi);
            auto const& vfrac_new =    EBFactory(lev).getVolFrac().const_array(mfi);

#ifdef AMREX_USE_MOVING_EB
	    //
	    // Redistribute
	    // (rho vel)^new = (rho vel)^old + dt * (
	    //                   div(rho vel u) + divtau_old + rho * f_t

	    // FIXME - I need to grow this box -- maybe safer to have 3 mfiters???
	    // incflo_compute_advective... does this same thing...
	    // except update now needs to be temporary, can't overwrite conv
	    Array4<Real      > const& vel_o = ld.velocity_o.array(mfi);
	    Array4<Real> const& update_v     = update.array(mfi);
	    Box const& gbx = amrex::grow(bx,ld.tracer.nGrow());
	    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		for (int n = 0; n < AMREX_SPACEDIM; ++n)
		{
		    vel_o(i,j,k,n) *= rho_old(i,j,k);
		}
	    });

	    // redistribute - own lambda...
	    redistribute_term(mfi, vel, update_v, vel_o, get_velocity_bcrec_device_ptr(), lev,
			      get_velocity_eb()[lev]->const_array(mfi));

	    // Don't divide out rho yet..
		
#else
	    
            if (m_diff_type == DiffusionType::Implicit) {

                if (use_tensor_correction)
                {
                    Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                    // Here divtau_o is the difference of tensor and scalar divtau_o!
                    if (m_advect_momentum) {
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                         vel(i,j,k,1) *= rho_old(i,j,k);,
                                         vel(i,j,k,2) *= rho_old(i,j,k););
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                        });
                    }
                } else {
                    if (m_advect_momentum) {
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                         vel(i,j,k,1) *= rho_old(i,j,k);,
                                         vel(i,j,k,2) *= rho_old(i,j,k););
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)););
                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)););
                        });
                    }
                }
            }
            else if (m_diff_type == DiffusionType::Crank_Nicolson)
            {

                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                if (m_advect_momentum) {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                     vel(i,j,k,1) *= rho_old(i,j,k);,
                                     vel(i,j,k,2) *= rho_old(i,j,k););
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0)+m_half*divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1)+m_half*divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)+m_half*divtau_o(i,j,k,2)););
                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+m_half*divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+m_half*divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+m_half*divtau_o(i,j,k,2)););
                    });
                }
            }
            else if (m_diff_type == DiffusionType::Explicit)
            {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                if (m_advect_momentum) {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        if (i==16 && j==4)//(vfrac_old(i,j,k) == 0.0 && vfrac_new(i,j,k)>0.0 ){
			{
                            Print()<<"vel pieces "<<vel(i,j,k,0)
                                <<" "<<rho_old(i,j,k)
                                <<" "<<dvdt(i,j,k,0)
                                <<" "<<rho_nph(i,j,k)
                                <<" "<<vel_f(i,j,k,0)
                                <<" "<<divtau_o(i,j,k,0)
                                <<" "<<rho_new(i,j,k)
                                <<std::endl;
                        }

                        AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                     vel(i,j,k,1) *= rho_old(i,j,k);,
                                     vel(i,j,k,2) *= rho_old(i,j,k););
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                    });
                }
            }
#endif
        } // mfi

#ifdef AMREX_USE_MOVING_EB
	//
	// Create S^m(rho*vel) and update-hat
	//
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	    Box const& bx = mfi.tilebox();

	    Array4<Real      > const& rhovel_o = ld.velocity_o.array(mfi);
	    Array4<Real      > const& rhovel   = ld.velocity.array(mfi);
	    // Store update-hat here...
	    Array4<Real      > const& dvdt_o   = ld.conv_velocity_o.array(mfi);

	    Box const& gbx = amrex::grow(bx,ld.velocity.nGrow());
	    FArrayBox update_fab(gbx,AMREX_SPACEDIM,The_Async_Arena());
	    Array4<Real> const& zero_arr = update_fab.array();
	    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		for (int n = 0; n < AMREX_SPACEDIM; ++n)
		{
		    zero_arr(i,j,k,n) = 0.;
		}
	    });

	    // redistribute - own lambda...
	    redistribute_term(mfi, dvdt_o, zero_arr, rhovel_o,
			      get_velocity_bcrec_device_ptr(), lev,
			      get_velocity_eb()[lev]->const_array(mfi));

	    // Make update-hat
	    // NOt saving S^m(rho*vel) for now, don't think we need it
	    // FIXME - I think this is okay as redistribute only uses vel on bx
	    Array4<Real const> const& rho      = ld.density.const_array(mfi);
	    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		for (int n = 0; n < AMREX_SPACEDIM; ++n)
		{
		    // Update-hat = S^m(rho*vel^n) - S^m(rho*vel^n - dt(A^n - D^n - f^n))
		    //            = S^m(rho*vel^n) - rhovel^(p,n+1)
		    dvdt_o(i,j,k,n) -= rhovel(i,j,k,n);
		    dvdt_o(i,j,k,n) /= l_dt;
			
		    // put back vel = rhovel/rho_new - non-ideal for MEB, but convection
		    // & proj expects just vel...
		    rhovel(i,j,k,n) /= rho(i,j,k); 
		}
	    });
	} // mfi
#endif
    } // lev

    // *************************************************************************************
    // Solve diffusion equation for u* but using eta_old at old time
    // *************************************************************************************
    if (m_diff_type == DiffusionType::Crank_Nicolson || m_diff_type == DiffusionType::Implicit)
    {
        const int ng_diffusion = 1;
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillphysbc_velocity(lev, new_time, m_leveldata[lev]->velocity, ng_diffusion);
            fillphysbc_density (lev, new_time, m_leveldata[lev]->density , ng_diffusion);
        }

        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_dt : m_half*m_dt;
        diffuse_velocity(get_velocity_new(), get_density_new(), GetVecOfConstPtrs(vel_eta), dt_diff);
    }

#ifdef INCFLO_USE_MOVING_EB
    // *************************************************************************************
    //
    // Update the moving geometry and arrays
    //
    // *************************************************************************************

    // Need to remake the LevelData for the NodalProjection which computes U^(n+1)
    RemakeWithNewGeometry();

    // *************************************************************************************
    // Make half-time density with updated EB
    // *************************************************************************************
    Vector<MultiFab> density_nph_neweb;
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        density_nph_neweb.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(), Factory(lev));
	MultiFab::Copy(density_nph_neweb[lev], density_nph_oldeb[lev], 0, 0, 1, 1);
    }
#endif

    // VisMF::Write(m_leveldata[0]->density,"dens");
    // VisMF::Write(density_nph_neweb[0],"rnph");
    // VisMF::Write(m_leveldata[0]->velocity,"vel");


    // WritePlotFile();
    // static int count = 0;
    // count++;
    //if (count > 0) Abort();

    // *************************************************************************************
    //
    // Project velocity field, update pressure
    //
    // *************************************************************************************
    ApplyProjection(GetVecOfConstPtrs(density_nph_neweb),
                    new_time,m_dt,incremental_projection);

    // static int count = 0;
    // count++;
    // if (count > 2) Abort();

}
