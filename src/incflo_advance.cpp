#include <incflo.H>

#ifdef AMREX_USE_MOVING_EB
#include <hydro_redistribution.H>
#endif

using namespace amrex;

void incflo::Advance()
{
    BL_PROFILE("incflo::Advance");

    // Start timing current time step
    Real strt_step = static_cast<Real>(ParallelDescriptor::second());

    // Compute time step size
    int initialisation = 0;
    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);
    ComputeDt(initialisation, explicit_diffusion);

    // Set new and old time to correctly use in fillpatching
    for(int lev = 0; lev <= finest_level; lev++)
    {
        m_t_old[lev] = m_cur_time;
        m_t_new[lev] = m_cur_time + m_dt;
    }

    if (m_verbose > 0)
    {
        amrex::Print() << "\nStep " << m_nstep + 1
                       << ": from old_time " << m_cur_time
                       << " to new time " << m_cur_time + m_dt
                       << " with dt = " << m_dt << ".\n" << std::endl;
    }

    // Note that fillpatch_xx(new_time) won't work to copy new data into old container.
    // That would be passing inconsistent info about the time of the old container
    copy_from_new_to_old_velocity();
    copy_from_new_to_old_density();
    copy_from_new_to_old_tracer();


    int ng = nghost_state();
    for (int lev = 0; lev <= finest_level; ++lev) {
        // Note: fillpatch pulls any EBFactory info from the coarse level
        fillpatch_velocity(lev, m_t_old[lev], m_leveldata[lev]->velocity_o, ng);
        fillpatch_density(lev, m_t_old[lev], m_leveldata[lev]->density_o, ng);
        if (m_advect_tracer) {
            fillpatch_tracer(lev, m_t_old[lev], m_leveldata[lev]->tracer_o, ng);
        }
    }


#ifdef INCFLO_USE_MOVING_EB
    // **********************************************************************************************
    //
    // Update the moving geometry and arrays
    //
    // **********************************************************************************************

//    VisMF::Write(m_leveldata[0]->density_o,"do1");

    // Create the time n+1 geometry and associated Factories
    MakeNewEBGeometry(m_t_new[0]);
    MakeFactoryWithNewGeometry();
#endif

    // FIXME - don;t know that we need this here. Shouldn't this be good from the
    // last time step? -- that's only true for moving EB probably
    // Also, need this to fill eb_vel for 1st time step with MEB
#ifdef AMREX_USE_EB
    if (m_eb_flow.enabled) {
       for (int lev = 0; lev <= finest_level; ++lev) {
         if (m_eb_flow.is_omega) {
            set_eb_velocity_for_rotation(lev, m_t_old[lev], *get_velocity_eb()[lev], 1);
         } else {
            set_eb_velocity(lev, m_t_old[lev], *get_velocity_eb()[lev], 1);
         }
         set_eb_density(lev, m_t_old[lev], *get_density_eb()[lev], 1);
         set_eb_tracer(lev, m_t_old[lev], *get_tracer_eb()[lev], 1);
       }
    }
#endif


#ifdef INCFLO_USE_MOVING_EB
    for (int lev = 0; lev <= finest_level; lev++)
    {
// This is needed for vel and tracer provided we continue with the same
        // strategy of passing back an update rather than a state...

        // FIXME - need some way to make sure target volfrac is consistent between
        // here and other calls
        Real target_volfrac = Redistribution::defaults::target_vol_fraction;
        // FOR varible density, we have to take the NU cell's merging neighbor
        // Needs one filled ghost cell. Fills only valid region
        Redistribution::FillNewlyUncovered(m_leveldata[lev]->velocity_o,
                                           OldEBFactory(lev), EBFactory(lev),
                                           *get_velocity_eb()[lev],
                                           geom[lev], target_volfrac);
        fillpatch_velocity(lev, m_t_old[lev], m_leveldata[lev]->velocity_o, ng);
        // Note: fillpatch pulls any EBFactory info from the coarse level, so
        // as long as EB stays contained within the finest level, we're fine,
        // but this probably won't work otherwise

        //FIXME
        // Update in ApplyPredictor assumes new vel is the same as old vel
        // Should think about whether to change ApplyPredictor, do this or
        // do a copy of vel_old
        Redistribution::FillNewlyUncovered(m_leveldata[lev]->velocity,
                                           OldEBFactory(lev), EBFactory(lev),
                                           *get_velocity_eb()[lev],
                                           geom[lev], target_volfrac);
        fillpatch_velocity(lev, m_t_old[lev], m_leveldata[lev]->velocity, ng);

        Redistribution::FillNewlyUncovered(m_leveldata[lev]->density_o,
                                           OldEBFactory(lev), EBFactory(lev),
                                           *get_velocity_eb()[lev],
                                           geom[lev], target_volfrac);
        fillpatch_density(lev, m_t_old[lev], m_leveldata[lev]->density_o, ng);

        if (m_ntrac > 0) {
            Redistribution::FillNewlyUncovered(m_leveldata[lev]->tracer_o,
                                               OldEBFactory(lev), EBFactory(lev),
                                               *get_velocity_eb()[lev],
                                               geom[lev], target_volfrac);
            fillpatch_tracer(lev, m_t_old[lev], m_leveldata[lev]->tracer_o, ng);
        }

        Redistribution::FillNewlyUncovered(m_leveldata[lev]->gp,
                                           OldEBFactory(lev), EBFactory(lev),
                                           *get_velocity_eb()[lev],
                                           geom[lev], target_volfrac);

// FIXME will we also need to worry about all the pieces of U*, forces, etc???
        // this should be done in pred/corr
    }
#endif


    ApplyPredictor();

        //FIXME
    // this will overwrite the previous time plotfile
    // WritePlotFile();
    // static int count=0; count++;
    // if (count>44) Abort();

    if (m_advection_type == "MOL") {
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillpatch_velocity(lev, m_t_new[lev], m_leveldata[lev]->velocity, ng);
            fillpatch_density(lev, m_t_new[lev], m_leveldata[lev]->density, ng);
            if (m_advect_tracer) {
                fillpatch_tracer(lev, m_t_new[lev], m_leveldata[lev]->tracer, ng);
            }
        }

    //     //FIXME
    // // this will overwrite the previous time plotfile
    // WritePlotFile();
    // static int count=0; count++;
    // //if (count>2) Abort();

        ApplyCorrector();
    }

#if 0
    // This sums over all levels
    if (m_test_tracer_conservation) {
        Real sum = volumeWeightedSum(get_tracer_new_const(),0,geom,ref_ratio);
        amrex::Print() << "Sum tracer volume wgt2 = " << m_cur_time+m_dt << " " <<
                           sum << std::endl;
    }
#endif

    // Stop timing current time step
    Real end_step = static_cast<Real>(ParallelDescriptor::second()) - strt_step;
    ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
    if (m_verbose > 0)
    {
        amrex::Print() << "Time per step " << end_step << std::endl;
    }
}

