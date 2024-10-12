#include <AMReX_BC_TYPES.H>
#include <incflo.H>
#ifdef AMREX_USE_EB
#include <AMReX_EB_Redistribution.H>
#endif

using namespace amrex;

void incflo::ReadParameters ()
{
    {
        // Variables without prefix in inputs file
    ParmParse pp;

    pp.query("stop_time", m_stop_time);
    pp.query("max_step", m_max_step);
    pp.query("steady_state", m_steady_state);
    }

    { // Prefix amr
        ParmParse pp("amr");
        pp.query("regrid_int", m_regrid_int);
#ifdef AMREX_USE_EB
        pp.query("refine_cutcells", m_refine_cutcells);
#endif
#ifdef INCFLO_USE_PARTICLES
        pp.query("refine_particles", m_refine_particles);
#endif
        pp.query("KE_int", m_KE_int);

    } // end prefix amr

    { // Prefix incflo
        ParmParse pp("incflo");

        pp.query("verbose", m_verbose);

        pp.query("steady_state_tol", m_steady_state_tol);
        pp.query("initial_iterations", m_initial_iterations);
        pp.query("do_initial_proj", m_do_initial_proj);
        pp.query("do_initial_pressure_proj", m_do_initial_pressure_proj);

        pp.query("fixed_dt", m_fixed_dt);
        pp.query("cfl", m_cfl);

        // This will multiply the time-step in the very first step only
        pp.query("init_shrink", m_init_shrink);
        if (m_init_shrink > 1.0) {
            amrex::Abort("We require m_init_shrink <= 1.0");
        }

        // This limits dt growth per time step
        pp.query("dt_change_max", m_dt_change_max);
        if ( m_dt_change_max < 1.0 || m_dt_change_max > 1.1 ) {
            amrex::Abort("We require 1. < dt_change_max <= 1.1");
        }

        // Physics
        pp.queryarr("delp", m_delp, 0, AMREX_SPACEDIM);
        pp.queryarr("gravity", m_gravity, 0, AMREX_SPACEDIM);

        pp.query("constant_density"         , m_constant_density);
        pp.query("advect_tracer"            , m_advect_tracer);
        pp.query("test_tracer_conservation" , m_test_tracer_conservation);

        // Are we advecting velocity or momentum (default is velocity)
        pp.query("advect_momentum"                  , m_advect_momentum);

        // Are we using MOL or Godunov?
        pp.query("advection_type"                   , m_advection_type);
        pp.query("use_ppm"                          , m_godunov_ppm);
        pp.query("godunov_use_forces_in_trans"      , m_godunov_use_forces_in_trans);
        pp.query("godunov_include_diff_in_forcing"  , m_godunov_include_diff_in_forcing);
        pp.query("use_mac_phi_in_godunov"           , m_use_mac_phi_in_godunov);
        pp.query("use_cc_proj"                      , m_use_cc_proj);

        // What type of redistribution algorithm;
        // {NoRedist, FluxRedist, StateRedist}
#ifdef AMREX_USE_EB
        pp.query("redistribution_type"              , m_redistribution_type);
        if (m_redistribution_type != "NoRedist" &&
            m_redistribution_type != "FluxRedist" &&
            m_redistribution_type != "StateRedist")
            amrex::Abort("redistribution type must be NoRedist, FluxRedist, or StateRedist");

        if (m_advection_type == "Godunov" && m_godunov_ppm) amrex::Abort("Can't use PPM with EBGodunov");
        pp.query("write_geom_chk", m_write_geom_chk);
#endif

        if (m_advection_type == "MOL") m_godunov_include_diff_in_forcing = false;

        if (m_advection_type != "MOL" && m_advection_type != "Godunov" && m_advection_type != "BDS")
            amrex::Abort("advection type must be MOL, Godunov, or BDS");

        // The default for diffusion_type is 2, i.e. the default m_diff_type is DiffusionType::Implicit
        int diffusion_type = 2;
        pp.query("diffusion_type", diffusion_type);
        if (diffusion_type == 0) {
            m_diff_type = DiffusionType::Explicit;
        } else if (diffusion_type == 1) {
            m_diff_type = DiffusionType::Crank_Nicolson;
        } else if (diffusion_type == 2) {
            m_diff_type = DiffusionType::Implicit;
        } else {
            amrex::Abort("We currently require diffusion_type = 0 for explicit, 1 for Crank-Nicolson or 2 for implicit");
        }

        // Default is true; should we use tensor solve instead of separate solves for each component?
        pp.query("use_tensor_solve",use_tensor_solve);
        pp.query("use_tensor_correction",use_tensor_correction);

        if (use_tensor_solve && use_tensor_correction) {
            amrex::Abort("We cannot have both use_tensor_solve and use_tensor_correction be true");
        }

        if (m_diff_type != DiffusionType::Implicit && use_tensor_correction) {
            amrex::Abort("We cannot have use_tensor_correction be true and diffusion type not Implicit");
        }

        if (m_advection_type == "MOL" && m_cfl > 0.5) {
            amrex::Abort("We currently require cfl <= 0.5 when using the MOL advection scheme");
        }
        if (m_advection_type != "MOL" && m_cfl > 1.0) {
            amrex::Abort("We currently require cfl <= 1.0 when using this advection scheme");
        }

        pp.query("ntrac", m_ntrac);

        if (m_ntrac <= 0) m_advect_tracer = false;

        if (m_ntrac < 1) {
            amrex::Abort("We currently require at least one tracer");
        }

        // Initial conditions
        pp.query("probtype", m_probtype);
        pp.query("ic_u", m_ic_u);
        pp.query("ic_v", m_ic_v);
        pp.query("ic_w", m_ic_w);
        pp.query("ic_p", m_ic_p);
        if ( !pp.queryarr("ic_t", m_ic_t, 0, m_ntrac) ) {
            m_ic_t.resize(m_ntrac, 0.);
        }

        // Viscosity (if constant)
        pp.query("mu", m_mu);

        // Density (if constant)
        pp.query("ro_0", m_ro_0);
        AMREX_ALWAYS_ASSERT(m_ro_0 >= 0.0);

        // Scalar diffusion coefficients
        m_mu_s.resize(m_ntrac, 0.0);
        pp.queryarr("mu_s", m_mu_s, 0, m_ntrac );

        amrex::Print() << "Scalar diffusion coefficients " << std::endl;
        for (int i = 0; i < m_ntrac; i++) {
            amrex::Print() << "Tracer diffusion coeff: " << i << ":" << m_mu_s[i] << std::endl;
        }
        //vof parameters
        pp.query("vof_advect_tracer", m_vof_advect_tracer);
        if (m_vof_advect_tracer){
           //the default of the density of VOF phase is same as the background fluid
           m_ro_s.resize(m_ntrac, m_ro_0);
           pp.queryarr("ro_s", m_ro_s, 0, m_ntrac );
           // the default of the surface tension is zero
           m_sigma.resize(m_ntrac, 0.);
           pp.queryarr("sigma", m_sigma, 0, m_ntrac );
        }
        if(m_vof_advect_tracer){
            m_update_density_from_vof = true;
            m_constant_density = false;
        }
        pp.query("number_of_averaging", m_number_of_averaging);

    } // end prefix incflo

    ReadIOParameters();
    ReadRheologyParameters();

    { // Prefix mac
        ParmParse pp_mac("mac_proj");
        pp_mac.query( "mg_rtol"                , m_mac_mg_rtol );
        pp_mac.query( "mg_atol"                , m_mac_mg_atol );
        pp_mac.query( "mg_max_coarsening_level", m_mac_mg_max_coarsening_level );
    } // end prefix mac

    { // Prefix nodal
        ParmParse pp_nodal("nodal_proj");
        pp_nodal.query( "mg_max_coarsening_level", m_nodal_mg_max_coarsening_level );
        pp_nodal.query( "mg_rtol"                , m_nodal_mg_rtol );
        pp_nodal.query( "mg_atol"                , m_nodal_mg_atol );
    } // end prefix nodal

#ifdef AMREX_USE_EB
    { // Prefix eb_flow
       ParmParse pp_eb_flow("eb_flow");

       pp_eb_flow.query("density", m_eb_flow.density);

       m_eb_flow.tracer.resize(m_ntrac, 0.0);
       pp_eb_flow.queryarr("tracer", m_eb_flow.tracer, 0, m_ntrac);

       if (pp_eb_flow.contains("vel_mag")) {
          m_eb_flow.enabled = true;
          m_eb_flow.is_mag = true;
          pp_eb_flow.query("vel_mag", m_eb_flow.vel_mag);
       } else if (pp_eb_flow.contains("velocity")) {
          m_eb_flow.enabled = true;
          pp_eb_flow.getarr("velocity", m_eb_flow.velocity, 0, AMREX_SPACEDIM);
       }

       if (pp_eb_flow.contains("normal")) {
          m_eb_flow.has_normal = true;
          pp_eb_flow.getarr("normal", m_eb_flow.normal, 0, AMREX_SPACEDIM);

          amrex::Real tol_deg(0.);
          pp_eb_flow.query("normal_tol", tol_deg);
          m_eb_flow.normal_tol = tol_deg*M_PI/amrex::Real(180.);
       }
    } // end prefix eb_flow
#endif

#ifdef INCFLO_USE_PARTICLES
    readTracerParticlesParams();
#endif

    if (m_use_cc_proj && max_level > 0) {
        amrex::Abort("Can't yet do multilevel with cell-centered projection");
    }
}

void incflo::ReadIOParameters()
{
    // Prefix amr
    ParmParse pp("amr");

    pp.query("check_file", m_check_file);
    pp.query("check_int", m_check_int);
    pp.query("restart", m_restart_file);

    pp.query("plotfile_on_restart", m_plotfile_on_restart);
    pp.query("regrid_on_restart", m_regrid_on_restart);

    pp.query("plot_file", m_plot_file);
    pp.query("plot_int"       , m_plot_int);
    pp.query("plot_per_exact" , m_plot_per_exact);
    pp.query("plot_per_approx", m_plot_per_approx);

    if ( (m_plot_int       > 0 && m_plot_per_exact  > 0) ||
         (m_plot_int       > 0 && m_plot_per_approx > 0) ||
         (m_plot_per_exact > 0 && m_plot_per_approx > 0) )
       amrex::Abort("Must choose only one of plot_int or plot_per_exact or plot_per_approx");

    // The plt_ccse_regtest resets the defaults,
    //     but we can over-ride those below
    int plt_ccse_regtest = 0;
    pp.query("plt_ccse_regtest", plt_ccse_regtest);

    if (plt_ccse_regtest != 0)
    {
        m_plt_velx       = 1;
        m_plt_vely       = 1;
        m_plt_velz       = 1;
        m_plt_gpx        = 1;
        m_plt_gpy        = 1;
        m_plt_gpz        = 1;
        m_plt_rho        = 1;
        m_plt_tracer     = 1;
        m_plt_p          = 0;
        m_plt_macphi     = 0;
        m_plt_eta        = 0;
        m_plt_vort       = 0;
        m_plt_magvel     = 0;
        m_plt_strainrate = 0;
        m_plt_divu       = 0;
        m_plt_vfrac      = 0;
#ifdef INCFLO_USE_PARTICLES
        m_plt_particle_count = 1;
#endif
    }

    // Which variables to write to plotfile

    pp.query("plt_velx",       m_plt_velx  );
    pp.query("plt_vely",       m_plt_vely  );
    pp.query("plt_velz",       m_plt_velz  );

    pp.query("plt_gpx",        m_plt_gpx );
    pp.query("plt_gpy",        m_plt_gpy );
    pp.query("plt_gpz",        m_plt_gpz );

    pp.query("plt_rho",        m_plt_rho   );
    pp.query("plt_tracer",     m_plt_tracer);
    pp.query("plt_p   ",       m_plt_p     );
    pp.query("plt_macphi",     m_plt_macphi);
    pp.query("plt_eta",        m_plt_eta   );
    pp.query("plt_magvel",     m_plt_magvel);
    pp.query("plt_vort",       m_plt_vort  );
    pp.query("plt_strainrate", m_plt_strainrate);
    pp.query("plt_divu",       m_plt_divu  );
    pp.query("plt_vfrac",      m_plt_vfrac );

    pp.query("plt_forcing",    m_plt_forcing );

    pp.query("plt_error_u",    m_plt_error_u );
    pp.query("plt_error_v",    m_plt_error_v );
    pp.query("plt_error_w",    m_plt_error_w );
    pp.query("plt_error_p",    m_plt_error_p );
    pp.query("plt_error_mac_p",m_plt_error_mac_p );

#ifdef INCFLO_USE_PARTICLES
    pp.query("plt_particle_count", m_plt_particle_count );
#endif
}

//
// Perform initial pressure iterations
//
void incflo::InitialIterations ()
{
    BL_PROFILE("incflo::InitialIterations()");

    copy_from_new_to_old_velocity();
    copy_from_new_to_old_density();
    copy_from_new_to_old_tracer();

    int initialisation = 1;
    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);
    ComputeDt(initialisation, explicit_diffusion);

    if (m_verbose && m_initial_iterations > 0)
    {
        amrex::Print() << "Doing initial pressure iterations with dt = " << m_dt << std::endl;
    }

    auto mac_phi = get_mac_phi();

    for (int lev = 0; lev <= finest_level; ++lev) m_t_old[lev] = m_t_new[lev];
    for (int lev = 0; lev <= finest_level; ++lev) mac_phi[lev]->setVal(0.);

    int ng = nghost_state();
    for (int lev = 0; lev <= finest_level; ++lev) {
            fillpatch_velocity(lev, m_t_old[lev], m_leveldata[lev]->velocity_o, ng);
        fillpatch_density(lev, m_t_old[lev], m_leveldata[lev]->density_o, ng);
        if (m_advect_tracer) {
            fillpatch_tracer(lev, m_t_old[lev], m_leveldata[lev]->tracer_o, ng);
        }
    }

    for (int iter = 0; iter < m_initial_iterations; ++iter)
    {
        if (m_verbose) amrex::Print() << "\n In initial_iterations: iter = " << iter << "\n";

     ApplyPredictor(true);

        copy_from_old_to_new_velocity();
        copy_from_old_to_new_density();
        copy_from_old_to_new_tracer();
    }

    // Reset dt to get initial step as specified, otherwise we can see increase to dt
    m_prev_dt = Real(-1.0);
    m_dt = Real(-1.0);
}

// Project velocity field to make sure initial velocity is divergence-free
void incflo::InitialProjection()
{
    BL_PROFILE("incflo::InitialProjection()");

    // *************************************************************************************
    // Allocate space for the temporary MAC velocities
    // *************************************************************************************
    Vector<MultiFab> u_mac_tmp(finest_level+1), v_mac_tmp(finest_level+1), w_mac_tmp(finest_level+1);
    int ngmac = nghost_mac();

    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(u_mac_tmp[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(0)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev));,
                     v_mac_tmp[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(1)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev));,
                     w_mac_tmp[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(2)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev)););
        if (ngmac > 0) {
            AMREX_D_TERM(u_mac_tmp[lev].setBndry(0.0);,
                         v_mac_tmp[lev].setBndry(0.0);,
                         w_mac_tmp[lev].setBndry(0.0););
        }
    }

    Real dummy_dt = 1.0;
    bool incremental_projection = false;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        m_leveldata[lev]->density.FillBoundary(geom[lev].periodicity());
    }

    ApplyProjection(get_density_new_const(),
                    AMREX_D_DECL(GetVecOfPtrs(u_mac_tmp), GetVecOfPtrs(v_mac_tmp),
                    GetVecOfPtrs(w_mac_tmp)),m_cur_time,dummy_dt,incremental_projection);


    // We set p and gp back to zero (p0 may still be still non-zero)
    for (int lev = 0; lev <= finest_level; lev++)
    {
        m_leveldata[lev]->p_nd.setVal(0.0);
        m_leveldata[lev]->p_cc.setVal(0.0);
        m_leveldata[lev]->gp.setVal(0.0);
    }
}

// Project to enforce hydrostatic equilibrium
void incflo::InitialPressureProjection()
{
    BL_PROFILE("incflo::InitialPressureProjection()");

    if (m_verbose > 0) { Print() << " Initial pressure projection \n"; }

    Real dummy_dt = 1.0;
    int  nGhost = 1;

    // fixme??? are density ghosts fill already???
    //I think we only need to worry about this if doing outflow bcs...
    for (int lev = 0; lev <= finest_level; lev++)
    {
        m_leveldata[lev]->density.FillBoundary(geom[lev].periodicity());
    }

    // *************************************************************************************
    // Allocate space for the temporary MAC velocities
    // *************************************************************************************
    Vector<MultiFab> u_mac_tmp(finest_level+1), v_mac_tmp(finest_level+1), w_mac_tmp(finest_level+1);
    int ngmac = nghost_mac();

    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(u_mac_tmp[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(0)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev));,
                     v_mac_tmp[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(1)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev));,
                     w_mac_tmp[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(2)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev)););
        if (ngmac > 0) {
            AMREX_D_TERM(u_mac_tmp[lev].setBndry(0.0);,
                         v_mac_tmp[lev].setBndry(0.0);,
                         w_mac_tmp[lev].setBndry(0.0););
        }
    }

    // Set the velocity to the gravity field
    Vector<MultiFab> vel(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, nGhost,
                        MFInfo(), *m_factory[lev]);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[lev].setVal(m_gravity[idim], idim, 1, 1);
        }

        auto& ld = *m_leveldata[lev];

        Real rho0 = m_ro_0;
        if (rho0 > 0.0) {
            for (MFIter mfi(ld.density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.growntilebox(1);
                Array4<Real const> const& rho_arr = ld.density.const_array(mfi);
                Array4<Real      > const& vel_arr = vel[lev].array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhofac = (rho_arr(i,j,k) - rho0) / rho_arr(i,j,k);
                    AMREX_D_TERM(vel_arr(i,j,k,0) *= rhofac;,
                                 vel_arr(i,j,k,1) *= rhofac;,
                                 vel_arr(i,j,k,2) *= rhofac;);
                });
            } // mfi
        } // rho0
    } // lev

    // Cell-centered divergence condition source term
    // Always zero this here
    Vector<MultiFab*> Source(finest_level+1, nullptr);

    // FIXME FIXME FIXME - THIS ONLY WORKS RIGHT FOR NODAL PROJ
    ApplyProjection(get_density_new_const(), GetVecOfPtrs(vel), Source,
                    m_cur_time, dummy_dt, false /*incremental*/,
                    true /*set_inflow_bc*/);
}

#ifdef AMREX_USE_EB
void
incflo::InitialRedistribution ()
{
    // Next we must redistribute the initial solution if we are going to use
    // StateRedist redistribution scheme
    if (m_redistribution_type == "StateRedist")
    {
      for (int lev = 0; lev <= finest_level; lev++)
      {
        auto& ld = *m_leveldata[lev];

        // We use the "old" data as the input here
        // We must fill internal ghost values before calling redistribution
        // We also need any physical boundary conditions imposed if we are
        //    calling state redistribution (because that calls the slope routine)

        ld.velocity.FillBoundary(geom[lev].periodicity());
        MultiFab::Copy(ld.velocity_o, ld.velocity, 0, 0, AMREX_SPACEDIM, ld.velocity.nGrow());
        fillpatch_velocity(lev, m_t_new[lev], ld.velocity_o, 3);

        if (!m_constant_density||m_vof_advect_tracer)
        {
            ld.density.FillBoundary(geom[lev].periodicity());
            MultiFab::Copy(ld.density_o, ld.density, 0, 0, 1, ld.density.nGrow());
            fillpatch_density(lev, m_t_new[lev], ld.density_o, 3);
        }
        if (m_advect_tracer)
        {
            ld.tracer.FillBoundary(geom[lev].periodicity());
            MultiFab::Copy(ld.tracer_o, ld.tracer, 0, 0, m_ntrac, ld.tracer.nGrow());
            fillpatch_tracer(lev, m_t_new[lev], ld.tracer_o, 3);
        }

        for (MFIter mfi(ld.density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& fact = EBFactory(lev);

            EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
            Array4<EBCellFlag const> const& flag = flagfab.const_array();

            if ( (flagfab.getType(bx)                != FabType::covered) &&
                 (flagfab.getType(amrex::grow(bx,4)) != FabType::regular) )
            {
                Array4<Real const> AMREX_D_DECL(fcx, fcy, fcz), ccc, vfrac, AMREX_D_DECL(apx, apy, apz);
                AMREX_D_TERM(fcx = fact.getFaceCent()[0]->const_array(mfi);,
                             fcy = fact.getFaceCent()[1]->const_array(mfi);,
                             fcz = fact.getFaceCent()[2]->const_array(mfi););
                ccc   = fact.getCentroid().const_array(mfi);
                AMREX_D_TERM(apx = fact.getAreaFrac()[0]->const_array(mfi);,
                             apy = fact.getAreaFrac()[1]->const_array(mfi);,
                             apz = fact.getAreaFrac()[2]->const_array(mfi););
                vfrac = fact.getVolFrac().const_array(mfi);

                int ncomp = AMREX_SPACEDIM;
                auto const& bc_vel = get_velocity_bcrec_device_ptr();
                ApplyInitialRedistribution( bx,ncomp,
                                          ld.velocity.array(mfi), ld.velocity_o.array(mfi),
                                          flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                          AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                          bc_vel, geom[lev], m_redistribution_type);

                if (!m_constant_density)
                {
                    ncomp = 1;
                    auto const& bc_den = get_density_bcrec_device_ptr();
                    ApplyInitialRedistribution( bx,ncomp,
                                              ld.density.array(mfi), ld.density_o.array(mfi),
                                              flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                              AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                              bc_den, geom[lev], m_redistribution_type);
                }
                if (m_advect_tracer)
                {
                    ncomp = m_ntrac;
                    auto const& bc_tra = get_tracer_bcrec_device_ptr();
                    ApplyInitialRedistribution( bx,ncomp,
                                              ld.tracer.array(mfi), ld.tracer_o.array(mfi),
                                              flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                              AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                              bc_tra, geom[lev], m_redistribution_type);
                }
            }
        }

        // We fill internal ghost values after calling redistribution
        ld.velocity.FillBoundary(geom[lev].periodicity());
        ld.density.FillBoundary(geom[lev].periodicity());
        ld.tracer.FillBoundary(geom[lev].periodicity());
    }
  }
}
#endif
