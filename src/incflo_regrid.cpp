#include <incflo.H>

using namespace amrex;

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void incflo::MakeNewLevelFromCoarse (int lev,
                                     Real time,
                                     const BoxArray& ba,
                                     const DistributionMapping& dm)
{
    BL_PROFILE("incflo::MakeNewLevelFromCoarse()");

    if (m_verbose > 0) {
        amrex::Print() << "Making new level " << lev << " from coarse" << std::endl;
    }

#ifdef AMREX_USE_EB
    std::unique_ptr<EBFArrayBoxFactory> new_fact = makeEBFabFactory(geom[lev], ba, dm,
                                                                    {nghost_eb_basic(),
                                                                    nghost_eb_volume(),
                                                                    nghost_eb_full()},
                                                                    EBSupport::full);
#else
    std::unique_ptr<FabFactory<FArrayBox> > new_fact(new FArrayBoxFactory());
#endif

    std::unique_ptr<LevelData> new_leveldata
        (new LevelData(ba, dm, *new_fact, m_ntrac, nghost_state(),
                       m_advection_type,
                       m_diff_type==DiffusionType::Implicit,
                       use_tensor_correction,
                       m_advect_tracer));

    fillcoarsepatch_velocity(lev, time, new_leveldata->velocity, 0);
    fillcoarsepatch_density(lev, time, new_leveldata->density, 0);
    if (m_ntrac > 0) {
        fillcoarsepatch_tracer(lev, time, new_leveldata->tracer, 0);
    }
    fillcoarsepatch_gradp(lev, time, new_leveldata->gp, 0);
    new_leveldata->p_nd.setVal(0.0);

    m_leveldata[lev] = std::move(new_leveldata);
#ifdef INCFLO_USE_MOVING_EB
    m_old_factory[lev] = std::move(m_new_factory[lev]);
    m_new_factory[lev] = std::move(new_fact);
#else
    m_factory[lev] = std::move(new_fact);
#endif
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void incflo::RemakeLevel (int lev, Real time, const BoxArray& ba,
                          const DistributionMapping& dm)
{
    BL_PROFILE("incflo::RemakeLevel()");

    if (m_verbose > 0) {
        amrex::Print() << "Remaking level " << lev << std::endl;
    }

#ifdef AMREX_USE_EB
    std::unique_ptr<EBFArrayBoxFactory> new_fact = makeEBFabFactory(geom[lev], ba, dm,
                                                                    {nghost_eb_basic(),
                                                                     nghost_eb_volume(),
                                                                     nghost_eb_full()},
                                                                     EBSupport::full);
#else
    std::unique_ptr<FabFactory<FArrayBox> > new_fact(new FArrayBoxFactory());
#endif
    std::unique_ptr<LevelData> new_leveldata
        (new LevelData(ba, dm, *new_fact, m_ntrac, nghost_state(),
                       m_advection_type,
                       m_diff_type==DiffusionType::Implicit,
                       use_tensor_correction,
                       m_advect_tracer));
    fillpatch_velocity(lev, time, new_leveldata->velocity, 0);
    fillpatch_density(lev, time, new_leveldata->density, 0);
    if (m_ntrac > 0) {
        fillpatch_tracer(lev, time, new_leveldata->tracer, 0);
    }
    fillpatch_gradp(lev, time, new_leveldata->gp, 0);
    new_leveldata->p_nd.setVal(0.0);

    m_leveldata[lev] = std::move(new_leveldata);
#ifdef INCFLO_USE_MOVING_EB
    m_old_factory[lev] = std::move(m_new_factory[lev]);
    m_new_factory[lev] = std::move(new_fact);
#else
    m_factory[lev] = std::move(new_fact);
#endif

    m_diffusion_tensor_op.reset();
    m_diffusion_scalar_op.reset();

#ifdef AMREX_USE_EB
    macproj.reset(new Hydro::MacProjector(Geom(0,finest_level),
                      MLMG::Location::FaceCentroid,  // Location of mac_vec
                      MLMG::Location::FaceCentroid,  // Location of beta
                      MLMG::Location::CellCenter  ) ); // Location of solution variable phi
#else
    macproj.reset(new Hydro::MacProjector(Geom(0,finest_level)));
#endif
}

#ifdef AMREX_USE_MOVING_EB
void incflo::MakeFactoryWithNewGeometry ()
{
    BL_PROFILE("incflo::MakeFactoryWithNewGeometry()");

    if (m_verbose > 0) {
        amrex::Print() << "Updating Factory with new geometry" << std::endl;
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
        m_new_factory[lev] = makeEBFabFactory(m_eb_new, geom[lev], grids[lev], dmap[lev],
                                              {nghost_eb_basic(),
                                               nghost_eb_volume(),
                                               nghost_eb_full()},
                                              EBSupport::full);
    }
}

// Remake an existing level with a new geometry but nothing else changed
void incflo::RemakeWithNewGeometry ()
{
    BL_PROFILE("incflo::RemakeWithNewGeometry()");

    if (m_verbose > 0) {
        amrex::Print() << "Remaking with new geometry" << std::endl;
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
        std::unique_ptr<LevelData> new_leveldata
            (new LevelData(grids[lev], dmap[lev], *m_new_factory[lev], m_ntrac, nghost_state(),
                           m_advection_type,
                           m_diff_type==DiffusionType::Implicit,
                           use_tensor_correction,
                           m_advect_tracer));

        // We call this from the predictor, so haven't put anything in these new time
        // containers yet
        MultiFab::Copy(new_leveldata->conv_velocity_o , m_leveldata[lev]->conv_velocity_o,
                       0,0,AMREX_SPACEDIM,m_leveldata[lev]->conv_velocity_o.nGrow());
        //MultiFab::Copy(new_leveldata->conv_velocity , m_leveldata[lev]->conv_velocity,0,0,AMREX_SPACEDIM,0);
        MultiFab::Copy(new_leveldata->conv_density_o , m_leveldata[lev]->conv_density_o,
                       0,0,1,m_leveldata[lev]->conv_density_o.nGrow());
        //MultiFab::Copy(new_leveldata->conv_density , m_leveldata[lev]->conv_density,0,0,1,0);

        // FIXME - the diffusive containers are only defined in certain circumstances...
        MultiFab::Copy(new_leveldata->divtau_o , m_leveldata[lev]->divtau_o,0,0,AMREX_SPACEDIM,0);
        //MultiFab::Copy(new_leveldata->divtau , m_leveldata[lev]->divtau,0,0,AMREX_SPACEDIM,0);
// laps is only defined if m_advect_tracer??
        //MultiFab::Copy(new_leveldata->laps_o , m_leveldata[lev]->laps_o,0,0,m_ntrac,0);
        //MultiFab::Copy(new_leveldata->laps , m_leveldata[lev]->laps,0,0,m_ntrac,0);

        // FIXME - Why FB when I could just copy the filled ghosts?
        // new_leveldata->conv_velocity_o.FillBoundary(geom[lev].periodicity());
        // new_leveldata->conv_density_o.FillBoundary(geom[lev].periodicity());
        // new_leveldata->divtau_o.FillBoundary(geom[lev].periodicity());
        // new_leveldata->laps_o.FillBoundary(geom[lev].periodicity());

        Real old_time = m_cur_time;
        Real new_time = m_cur_time + m_dt;


// FIXME? Why FP here when we could just copy? Are ghost cells good? Are they needed?

        // FIXME? fillpatch uses m_leveldata (potentially old EB) under the covers,
        // and eb_cell_cons interpolator which pulls EB info from the coarse level
        // If we always call on level 0 first and go to finer, then it's probably okay...
        fillpatch_velocity(lev, old_time, new_leveldata->velocity_o, nghost_state());
        fillpatch_velocity(lev, new_time, new_leveldata->velocity, nghost_state());
        fillpatch_density(lev, old_time, new_leveldata->density_o, nghost_state());
        fillpatch_density(lev, new_time, new_leveldata->density, nghost_state());
        if (m_ntrac > 0) {
            fillpatch_tracer(lev, old_time, new_leveldata->tracer_o, nghost_state());
            fillpatch_tracer(lev, new_time, new_leveldata->tracer, nghost_state());
        }


        // time is really a dummy variable here, since we only carry one gradp (at
        // time n-1/2).
        fillpatch_gradp(lev, old_time, new_leveldata->gp, 0);
        new_leveldata->p_nd.setVal(0.0);

#if 0
        // Don't do this. Newly covered cell gets into calculation...
        // Let's fill the newly covered cells with 1e45 to be different
        EB_set_covered( new_leveldata->velocity  , 1.e45);
        EB_set_covered( new_leveldata->velocity_o, 1.e45);
        EB_set_covered( new_leveldata->density   , 1.e45);
        EB_set_covered( new_leveldata->density_o , 1.e45);
        if (m_ntrac > 0 && m_advect_tracer) {
            EB_set_covered( new_leveldata->tracer    , 1.e45);
            EB_set_covered( new_leveldata->tracer_o  , 1.e45);
        }
        EB_set_covered( new_leveldata->gp , 1.e45);

        EB_set_covered( new_leveldata->conv_velocity_o, 1.e45);
#endif

        // Update the member variable with the newly created version
        m_leveldata[lev] = std::move(new_leveldata);

#ifdef AMREX_USE_MOVING_EB
        // FIXME? Is this the right place to update eb_vals?
        // Is it better to recompute old or copy?
        if (m_eb_flow.enabled)
        {
            if (m_verbose >0) { Print()<<"Updating the eb_velocity..."<<std::endl; }

            if (m_eb_flow.is_omega) {
                set_eb_velocity_for_rotation(lev, m_t_old[lev], *get_velocity_eb(old_time)[lev],
                                             get_velocity_eb(old_time)[lev]->nGrow());
                set_eb_velocity_for_rotation(lev, m_t_new[lev], *get_velocity_eb(new_time)[lev],
                                             get_velocity_eb(new_time)[lev]->nGrow());
            } else {
                set_eb_velocity(lev, m_t_old[lev], *get_velocity_eb(old_time)[lev],
                                get_velocity_eb(old_time)[lev]->nGrow());
                set_eb_velocity(lev, m_t_new[lev], *get_velocity_eb(new_time)[lev],
                                get_velocity_eb(new_time)[lev]->nGrow());
            }
            set_eb_density(lev, m_t_old[lev], *get_density_eb()[lev],
                           get_density_eb()[lev]->nGrow());
            set_eb_density(lev, m_t_new[lev], *get_density_eb()[lev],
                           get_density_eb()[lev]->nGrow());
            set_eb_tracer(lev, m_t_old[lev], *get_tracer_eb()[lev],
                          get_tracer_eb()[lev]->nGrow());
            set_eb_tracer(lev, m_t_new[lev], *get_tracer_eb()[lev],
                          get_tracer_eb()[lev]->nGrow());
        }
#endif

        // MATT -- reset macproj
        // FIXME - do we really want to do this here? Not sure the most logical place...
        macproj.reset(new Hydro::MacProjector(Geom(0,finest_level),
                                              MLMG::Location::FaceCentroid,  // Location of mac_vec
                                              MLMG::Location::FaceCentroid,  // Location of beta
                                              MLMG::Location::CellCenter  ) ); // Location of solution variable phi
    }
}
#endif

// Delete level data
// overrides the pure virtual function in AmrCore
void incflo::ClearLevel (int lev)
{
    BL_PROFILE("incflo::ClearLevel()");
    m_leveldata[lev].reset();
#ifdef INCFLO_USE_MOVING_EB
    m_new_factory[lev].reset();
    m_old_factory[lev].reset();
#else
    m_factory[lev].reset();
#endif
    m_diffusion_tensor_op.reset();
    m_diffusion_scalar_op.reset();
    macproj.reset();
}
