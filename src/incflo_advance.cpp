#include <incflo.H>

using namespace amrex;

void incflo::Advance()
{
    BL_PROFILE("incflo::Advance");

    // Start timing current time step
    Real start_step = static_cast<Real>(ParallelDescriptor::second());

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

    // **********************************************************************************************
    //
    // CRYO-PLUNGING: set velocity to be zero
    //
    // **********************************************************************************************
    // TODO: figure out how many hard resets are needed
    if (m_advect_energy || m_advect_tracer)
    {
        cryo_set_zero_vel();
        cryo_update_thermal();
    }
    
    copy_from_new_to_old_velocity();
    copy_from_new_to_old_density();
    copy_from_new_to_old_tracer();
    if (m_advect_energy) copy_from_new_to_old_energy();


    int ng = nghost_state();
    for (int lev = 0; lev <= finest_level; ++lev) {
        fillpatch_velocity(lev, m_t_old[lev], m_leveldata[lev]->velocity_o, ng);
        fillpatch_density(lev, m_t_old[lev], m_leveldata[lev]->density_o, ng);
        if (m_advect_tracer) {
            fillpatch_tracer(lev, m_t_old[lev], m_leveldata[lev]->tracer_o, ng);
        }
        if (m_advect_energy) {
            fillpatch_density(lev, m_t_old[lev], m_leveldata[lev]->one, ng);
            fillpatch_energy(lev, m_t_old[lev], m_leveldata[lev]->energy_o, ng);
            fillpatch_temp(lev, m_t_old[lev], m_leveldata[lev]->temp_o, ng);
        }
    }

#ifdef AMREX_USE_EB
    if (m_eb_flow.enabled) {
       for (int lev = 0; lev <= finest_level; ++lev) {
         set_eb_velocity(lev, m_t_old[lev], *get_velocity_eb()[lev], 1);
         set_eb_density(lev, m_t_old[lev], *get_density_eb()[lev], 1);
         set_eb_tracer(lev, m_t_old[lev], *get_tracer_eb()[lev], 1);
       }
    }
#endif

    ApplyPredictor();

    if (m_advection_type == "MOL") {
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillpatch_velocity(lev, m_t_new[lev], m_leveldata[lev]->velocity, ng);
            fillpatch_density(lev, m_t_new[lev], m_leveldata[lev]->density, ng);
            if (m_advect_tracer) {
                fillpatch_tracer(lev, m_t_new[lev], m_leveldata[lev]->tracer, ng);
            }
        }

        ApplyCorrector();
    }

    if (m_advect_energy || m_advect_tracer)
    {
        cryo_set_zero_vel();
    }

#if 0
    // This sums over all levels
    if (m_test_tracer_conservation) {
        Real sum = volumeWeightedSum(get_tracer_new_const(),0,geom,ref_ratio);
        amrex::Print() << "Sum tracer volume wgt2 = " << m_cur_time+m_dt << " " <<
                           sum << std::endl;
    }
#endif

    // **********************************************************************************************
    //
    // CRYO-PLUNGING: convect energy field
    //
    // **********************************************************************************************
    if (m_advect_energy)
    {
        // cryo_set_zero_vel(); // TODO: Check how many hard resets are needed
        // cryo_update_thermal();
        // cryo_compute_temp();
        // cryo_convect_energy();
        cryo_compute_temp();
        // HACK: Write energy and temperature directly
        cryo_copy_from_temp_to_tracer();
    }

    // Stop timing current time step
    Real end_step = static_cast<Real>(ParallelDescriptor::second()) - start_step;
    ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
    if (m_verbose > 0)
    {
        amrex::Print() << "Time per step " << end_step << std::endl;
    }
}

