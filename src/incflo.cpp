#include <incflo.H>

// Need this for TagCutCells
#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#include <utility>
#endif

#include <memory>

using namespace amrex;

incflo::incflo ()
{
    // NOTE: Geometry on all levels has just been defined in the AmrCore
    // constructor. No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    // Read inputs file using ParmParse
    ReadParameters();

#ifdef AMREX_USE_EB
    // This is needed before initializing level MultiFab
    MakeEBGeometry();
#endif

    // Initialize memory for data-array internals
    ResizeArrays();

    init_bcs();

    init_advection();

    set_background_pressure();
}

incflo::~incflo ()
= default;

void incflo::InitData ()
{
    BL_PROFILE("incflo::InitData()");

    if (m_restart_file.empty())
    {
        // This tells the AmrMesh class not to iterate when creating the initial
        // grid hierarchy
        // SetIterateToFalse();

        // This tells the Cluster routine to use the new chopping routine which
        // rejects cuts if they don't improve the efficiency
        SetUseNewChop();

        // This is an AmrCore member function which recursively makes new levels
        // with MakeNewLevelFromScratch.
        InitFromScratch(m_cur_time);

#ifdef AMREX_USE_EB
#ifdef INCFLO_USE_PARTICLES
        const auto& ebfact = EBFactory(0);
#endif
#endif

#ifdef INCFLO_USE_PARTICLES
        initializeTracerParticles( (ParGDBBase*)GetParGDB()
#ifdef AMREX_USE_EB
                                  ,ebfact
#endif
                                 );
#endif

#ifdef AMREX_USE_EB
        InitialRedistribution();
#endif

        if (m_do_initial_proj) {
            InitialProjection();

            if (m_do_initial_pressure_proj) {
                InitialPressureProjection();
            }
        }

        InitialIterations();

        // Set m_nstep to 0 before entering time loop
        m_nstep = 0;

        // xxxxx TODO averagedown ???

        if (m_check_int > 0) { WriteCheckPointFile(); }

        // Plot initial distribution
        if (m_plot_int > 0 || m_plot_per_exact > 0 || m_plot_per_approx > 0)
        {
            WritePlotFile();
            m_last_plt = 0;
        }
        if (m_KE_int > 0)
        {
            amrex::Abort("xxxxx m_KE_int todo");
//          amrex::Print() << "Time, Kinetic Energy: " << m_cur_time << ", " << ComputeKineticEnergy() << std::endl;
        }
    }
    else
    {
        // Read starting configuration from chk file.
        ReadCheckpointFile();

#ifdef INCFLO_USE_PARTICLES
        particleData.Redistribute();
#endif

        if (m_plotfile_on_restart)
        {
            WritePlotFile();
            m_last_plt = 0;
        }
    }

#ifdef AMREX_USE_EB
#if (AMREX_SPACEDIM == 3)
    ParmParse pp("incflo");
    bool write_eb_surface = false;
    pp.query("write_eb_surface", write_eb_surface);
    if (write_eb_surface) WriteMyEBSurface();
#endif
#endif

    if (m_verbose > 0 && ParallelDescriptor::IOProcessor()) {
        printGridSummary(amrex::OutStream(), 0, finest_level);
    }
}

void incflo::Evolve()
{
    BL_PROFILE("incflo::Evolve()");

    bool do_not_evolve = ((m_max_step == 0) || ((m_stop_time >= 0.) && (m_cur_time > m_stop_time)) ||
                           ((m_stop_time <= 0.) && (m_max_step <= 0)) || (m_max_step >= 0 && m_nstep >= m_max_step) )
                         && !m_steady_state;

    while(!do_not_evolve)
    {
        if (m_verbose > 0)
        {
            amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";
        }

        if (m_regrid_int > 0 && m_nstep > 0 && m_nstep%m_regrid_int == 0)
        {
            if (m_verbose > 0) amrex::Print() << "Regridding...\n";
            regrid(0, m_cur_time);
            if (m_verbose > 0 && ParallelDescriptor::IOProcessor()) {
                printGridSummary(amrex::OutStream(), 0, finest_level);
            }
        }

        // Advance to time t + dt
        Advance();
        m_nstep++;
        m_cur_time += m_dt;

        if (writeNow())
        {
            WritePlotFile();
            m_last_plt = m_nstep;
        }

        if(m_check_int > 0 && (m_nstep % m_check_int == 0))
        {
            WriteCheckPointFile();
            m_last_chk = m_nstep;
        }

        if(m_KE_int > 0 && (m_nstep % m_KE_int == 0))
        {
            amrex::Print() << "Time, Kinetic Energy: " << m_cur_time << ", " << ComputeKineticEnergy() << std::endl;
        }

        // Mechanism to terminate incflo normally.
        do_not_evolve = (m_steady_state && SteadyStateReached()) ||
                        ((m_stop_time > 0. && (m_cur_time >= m_stop_time - 1.e-12 * m_dt)) ||
                         (m_max_step >= 0 && m_nstep >= m_max_step));
    }

    // Output at the final time
    if( m_check_int > 0 && m_nstep != m_last_chk) {
        WriteCheckPointFile();
    }
    if( (m_plot_int > 0 || m_plot_per_exact > 0 || m_plot_per_approx > 0)
        && m_nstep != m_last_plt)
    {
        WritePlotFile();
    }
}

void
incflo::ApplyProjection (Vector<MultiFab const*> density,
                         Real time, Real scaling_factor, bool incremental)
{
    BL_PROFILE("incflo::ApplyProjection");
    ApplyNodalProjection(std::move(density),time,scaling_factor,incremental);
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void incflo::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& new_grids,
                                      const DistributionMapping& new_dmap)
{
    BL_PROFILE("incflo::MakeNewLevelFromScratch()");

    if (m_verbose > 0)
    {
        amrex::Print() << "Making new level " << lev << " from scratch" << std::endl;
        if (m_verbose > 2) {
            amrex::Print() << "with BoxArray " << new_grids << std::endl;
        }
    }

    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

#ifdef AMREX_USE_EB
    m_factory[lev] = makeEBFabFactory(geom[lev], grids[lev], dmap[lev],
                                      {nghost_eb_basic(),
                                       nghost_eb_volume(),
                                       nghost_eb_full()},
                                       EBSupport::full);
#else
    m_factory[lev] = std::make_unique<FArrayBoxFactory>();
#endif

    m_leveldata[lev] = std::make_unique<LevelData>(grids[lev], dmap[lev], *m_factory[lev],
                                                   this);

    m_t_new[lev] = time;
    m_t_old[lev] = time - Real(1.e200);

    if (m_restart_file.empty()) {
        prob_init_fluid(lev);
    }
    //make_mixedBC_mask(lev, grids[lev], dmap[lev]);

#ifdef AMREX_USE_EB
    macproj = std::make_unique<Hydro::MacProjector>(Geom(0,finest_level),
                      MLMG::Location::FaceCentroid,  // Location of mac_vec
                      MLMG::Location::FaceCentroid,  // Location of beta
                      MLMG::Location::CellCenter  ); // Location of solution variable phi
#else
    macproj = std::make_unique<Hydro::MacProjector>(Geom(0,finest_level));
#endif
}

bool
incflo::writeNow()
{
    bool write_now = false;

    if ( ( m_plot_int > 0 && (m_nstep % m_plot_int == 0) ) ||
         (m_plot_per_exact  > 0 && (std::abs(std::remainder(m_cur_time, m_plot_per_exact)) < 1.e-12) ) ) {
        write_now = true;
    } else if (m_plot_per_approx > 0.0) {
        // Check to see if we've crossed a m_plot_per_approx interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = static_cast<int>(std::round((m_cur_time-m_dt) / m_plot_per_approx));
        int num_per_new = static_cast<int>(std::round((m_cur_time     ) / m_plot_per_approx));

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next m_plot_per_approx interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * Real(10.0) * std::abs(m_cur_time);
        const Real next_plot_time = (num_per_old + 1) * m_plot_per_approx;

        if ((num_per_new == num_per_old) && std::abs(m_cur_time - next_plot_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && std::abs((m_cur_time - m_dt) - next_plot_time) <= eps)
            num_per_old += 1;

        if (num_per_old != num_per_new)
            write_now = true;
    }

    return write_now;
}
