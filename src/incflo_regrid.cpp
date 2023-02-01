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
    std::unique_ptr<EBFabFactory<FArrayBox> > new_fact(new FArrayBoxFactory());
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

#ifdef AMREX_USE_EB
#ifdef INCFLO_USE_MOVING_EB
void incflo::MakeNewGeometry (int lev, Real time)
{
    // Erase old EB
    EB2::IndexSpace::erase(const_cast<EB2::IndexSpace*>(m_eb_old));

    // Build a new EB
    MakeEBGeometry(time);

    m_eb_old = m_eb_new;
    m_eb_new = &(EB2::IndexSpace::top());

    m_old_factory[lev] = std::move(m_new_factory[lev]);

    m_new_factory[lev] = makeEBFabFactory(m_eb_new, geom[lev], grids[lev], dmap[lev],
                                          {nghost_eb_basic(),
                                           nghost_eb_volume(),
                                           nghost_eb_full()},
                                          EBSupport::full);

    EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(m_eb_new));
}

// Remake an existing level with a new geometry but nothing else changed
//FIXME - time parameter is misleading, remove...
void incflo::RemakeLevelWithNewGeometry (int lev, Real time)
{
    BL_PROFILE("incflo::RemakeLevelWithNewGeometry()");

    if (m_verbose > 0) {
        amrex::Print() << "Remaking level " << lev << " with new geometry" << std::endl;
    }

    // This has been called at the beginning of the timestep
    // MakeNewGeometry(lev,time);

    std::unique_ptr<LevelData> new_leveldata
        (new LevelData(grids[lev], dmap[lev], *m_new_factory[lev], m_ntrac, nghost_state(),
                       m_advection_type,
                       m_diff_type==DiffusionType::Implicit,
                       use_tensor_correction,
                       m_advect_tracer));

    // MultiFab::Copy(new_leveldata->velocity  , m_leveldata[lev]->velocity  ,0,0,AMREX_SPACEDIM,0);
    // MultiFab::Copy(new_leveldata->velocity_o, m_leveldata[lev]->velocity_o,0,0,AMREX_SPACEDIM,0);
    // MultiFab::Copy(new_leveldata->density   , m_leveldata[lev]->density  ,0,0,1,0);
    // MultiFab::Copy(new_leveldata->density_o , m_leveldata[lev]->density_o,0,0,1,0);
    // if (m_ntrac > 0) {
    //     MultiFab::Copy(new_leveldata->tracer  , m_leveldata[lev]->tracer  ,0,0,1,0);
    //     MultiFab::Copy(new_leveldata->tracer_o, m_leveldata[lev]->tracer_o,0,0,1,0);
    // }
    // MultiFab::Copy(new_leveldata->gp   , m_leveldata[lev]->gp   ,0,0,AMREX_SPACEDIM,0);
    // MultiFab::Copy(new_leveldata->p_nd , m_leveldata[lev]->p_nd ,0,0,1,0);

    MultiFab::Copy(new_leveldata->conv_velocity_o , m_leveldata[lev]->conv_velocity_o,0,0,AMREX_SPACEDIM,0);
    MultiFab::Copy(new_leveldata->conv_density_o , m_leveldata[lev]->conv_density_o,0,0,1,0);
    MultiFab::Copy(new_leveldata->conv_velocity , m_leveldata[lev]->conv_velocity,0,0,AMREX_SPACEDIM,0);
    MultiFab::Copy(new_leveldata->conv_density , m_leveldata[lev]->conv_density,0,0,1,0);

    VisMF::Write(m_leveldata[0]->density_o,"do10");
    // Fill in ghost cells for new MultiFabs (Matt - Not sure if 4 is the correct ng)
// This doesn;t work as expected because it relies on m_leveldata, which still has the old EB!!
    // however, the periodic fill is just supposed to look at whatever is in the valid region...
    // not sure why this doesn't work...
    Real old_time = m_cur_time;
    Real new_time = m_cur_time + m_dt;
    
    fillpatch_velocity(lev, old_time, new_leveldata->velocity_o, nghost_state());
    fillpatch_velocity(lev, new_time, new_leveldata->velocity, nghost_state());
    fillpatch_density(lev, old_time, new_leveldata->density_o, nghost_state());
    fillpatch_density(lev, new_time, new_leveldata->density, nghost_state());
    if (m_ntrac > 0) {
        fillpatch_tracer(lev, old_time, new_leveldata->tracer_o, nghost_state());
        fillpatch_tracer(lev, new_time, new_leveldata->tracer, nghost_state());
    }
    // time is really a dummy variable here. Since we only have one gp, FP will just take that
    fillpatch_gradp(lev, time, new_leveldata->gp, 0);
    new_leveldata->p_nd.setVal(0.0);
    VisMF::Write(new_leveldata->density_o,"do11");
    // No, we need to retain vals in NU cells in both new and old
#if 1
    // This should be okay to do...
    // Let's fill the newly covered cells with 1e45 to be different
    EB_set_covered( new_leveldata->velocity  , 1.e45);
    EB_set_covered( new_leveldata->velocity_o, 1.e45);
    EB_set_covered( new_leveldata->density   , 1.e45);
    EB_set_covered( new_leveldata->density_o , 1.e45);
    if (m_ntrac > 0) {
        EB_set_covered( new_leveldata->tracer    , 1.e45);
        EB_set_covered( new_leveldata->tracer_o  , 1.e45);
    }
    EB_set_covered( new_leveldata->gp , 1.e45);

    EB_set_covered( new_leveldata->conv_velocity_o, 1.e45);
#endif

    VisMF::Write(new_leveldata->density,"do12");
    
#if 0
    //FIXME - is this what we want to do for MOL pred-corr. new has been filled already from
    // update with MSRD, and we wouldn't want to re-define what rho_old is...
    //
    // Now let's make sure to fill cells that were previously covered but are now cut cell
    //amrex::Print() << "Fill Velocity" << std::endl;
    EB_fill_uncovered(lev,new_leveldata->velocity  , m_leveldata[lev]->velocity  );
    EB_fill_uncovered(lev,new_leveldata->velocity_o, m_leveldata[lev]->velocity_o);

    //amrex::Print() << "\nFill density" << std::endl;
    EB_fill_uncovered(lev,new_leveldata->density   , m_leveldata[lev]->density   );
    EB_fill_uncovered(lev,new_leveldata->density_o , m_leveldata[lev]->density_o );

    if (m_ntrac > 0) {
        EB_fill_uncovered(lev,new_leveldata->tracer    , m_leveldata[lev]->tracer    );
        EB_fill_uncovered(lev,new_leveldata->tracer_o  , m_leveldata[lev]->tracer_o  );
    }

    //amrex::Print() << "\nFill gp" << std::endl;
    EB_fill_uncovered(lev,new_leveldata->gp      , m_leveldata[lev]->gp      );

    //amrex::Print() << "\nFill p_nd" << std::endl;
    EB_fill_uncovered(lev,new_leveldata->p_nd    , m_leveldata[lev]->p_nd    );
#endif

    m_leveldata[lev] = std::move(new_leveldata);
    
    // MATT -- reset macproj
    macproj.reset(new Hydro::MacProjector(Geom(0,finest_level),
                      MLMG::Location::FaceCentroid,  // Location of mac_vec
                      MLMG::Location::FaceCentroid,  // Location of beta
                      MLMG::Location::CellCenter  ) ); // Location of solution variable phi


}
#endif
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

void incflo::EB_fill_uncovered (int lev, MultiFab& mf_new, MultiFab& mf_old)
{
    auto const& vfrac_old = OldEBFactory(lev).getVolFrac();
    auto const& vfrac_new =    EBFactory(lev).getVolFrac();

    for (MFIter mfi(mf_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real>       const& fab_new = mf_new.array(mfi);
        Array4<Real const> const& fab_old = mf_old.const_array(mfi);

        Array4<Real const> const&  vf_old = vfrac_old.const_array(mfi);
        Array4<Real const> const&  vf_new = vfrac_new.const_array(mfi);

        const int ncomp = mf_new.nComp();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // If new cell is uncovered... avg from neighbors that are cut or regular
            if (vf_old(i,j,k) == 0.0 && vf_new(i,j,k) > 0.0)
            {
                //amrex::Print() << "Need to fill cell " << IntVect(AMREX_D_DECL(i,j,k)) << std::endl;
                for (int n = 0; n < ncomp; n++)
                {
                    fab_new(i,j,k,n) = 0.;
                    Real den = 0.;

                    if (vf_old(i+1,j,k) > 0.0)
                    {
                        fab_new(i,j,k,n) += fab_old(i+1,j,k,n);
                        //amrex::Print() << "right fill: " << fab_old(i+1,j,k,n) << std::endl;
                        den += 1.;
                    }
                    if (vf_old(i-1,j,k) > 0.0)
                    {
                        fab_new(i,j,k,n) += fab_old(i-1,j,k,n);
                        //amrex::Print() << "left fill: " << fab_old(i-1,j,k,n) << std::endl;
                        den += 1.;
                    }
                    if (vf_old(i,j+1,k) > 0.0)
                    {
                        fab_new(i,j,k,n) += fab_old(i,j+1,k,n);
                        //amrex::Print() << "top fill: " << fab_old(i,j+1,k,n) << std::endl;
                        den += 1.;
                    }
                    if (vf_old(i,j-1,k) > 0.0)
                    {
                        fab_new(i,j,k,n) += fab_old(i,j-1,k,n);
                        //amrex::Print() << "bottom fill: " << fab_old(i,j-1,k,n) << std::endl;
                        den += 1.;
                    }
#if (AMREX_SPACEDIM == 3)
                    if (vf_old(i,j,k+1) > 0.0)
                    {
                        fab_new(i,j,k,n) += fab_old(i,j,k+1,n);
                        //amrex::Print() << "up fill: " << fab_old(i,j,k+1,n) << std::endl;
                        den += 1.;
                    }
                    if (vf_old(i,j,k-1) > 0.0)
                    {
                        fab_new(i,j,k,n) += fab_old(i,j,k-1,n);
                        //amrex::Print() << "down fill: " << fab_old(i,j,k-1,n) << std::endl;
                        den += 1.;
                    }
#endif

                    fab_new(i,j,k,n) = fab_new(i,j,k,n) / den;
                }
            }
        });
    }
}

void incflo::EB_fill_uncovered_with_zero (int lev, MultiFab& mf)
{
    auto const& vfrac_old = OldEBFactory(lev).getVolFrac();
    auto const& vfrac_new =    EBFactory(lev).getVolFrac();

    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real>       const& fab = mf.array(mfi);

        Array4<Real const> const&  vf_old = vfrac_old.const_array(mfi);
        Array4<Real const> const&  vf_new = vfrac_new.const_array(mfi);

        const int ncomp = mf.nComp();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // If new cell is uncovered... avg from neighbors that are cut or regular
            if (vf_old(i,j,k) == 0.0 && vf_new(i,j,k) > 0.0)
            {
                //amrex::Print() << "Need to fill cell with zero " << IntVect(AMREX_D_DECL(i,j,k)) << std::endl;
                for (int n = 0; n < ncomp; n++)
                {
		    // FIXME- for now make this not identically zero so inv does not cause error
                    fab(i,j,k,n) = 1.e-40;
                }
            }
        });
    }
}
