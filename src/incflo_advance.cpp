#include <incflo.H>
#include <volWgtSum.H>

using namespace amrex;

void incflo::Advance(Real orig_mass, Real& prev_mass)
{
    BL_PROFILE("incflo::Advance");

    // Start timing current time step
    Real strt_step = static_cast<Real>(ParallelDescriptor::second());

#ifdef INCFLO_USE_MOVING_EB
    // Point to the correct EB & Factory for current time
    EB2::IndexSpace::erase(const_cast<EB2::IndexSpace*>(m_eb_old));
    m_eb_old = m_eb_new;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        m_old_factory[lev] = std::move(m_new_factory[lev]);
    }
#endif


    // Compute time step size
    int initialisation = 0;
    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);
    ComputeDt(initialisation, explicit_diffusion);

    // Set new and old time to correctly use in fillpatching
    for(int lev = 0; lev <= finest_level; lev++)
    {
        m_t_old[lev] = m_cur_time;
        m_t_new[lev] = m_cur_time + m_dt;
	Print()<<"ADVANCE times : "<<m_t_old[lev]<<" "<<m_t_new[lev]<<std::endl;
    }

#ifdef INCFLO_USE_MOVING_EB
    // **********************************************************************************************
    //
    // Update the moving geometry and arrays
    //
    // **********************************************************************************************

    // Create the time n+1 geometry and associated Factories.
    // This moves m_old_factory to point to the correct EB.
    MakeNewEBGeometry(m_t_new[0]);
    MakeFactoryWithNewGeometry();
#endif

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


    ApplyPredictor();

#if 1
    int my_lev = 0;
#ifdef AMREX_USE_EB
    auto const& fact = EBFactory(my_lev);
    Real sum = volWgtSum(my_lev,get_density_new_const()[my_lev],0,fact);
#else
    Real sum = volWgtSum(my_lev,get_density_new_const()[my_lev],0);
#endif

    auto const dx = geom[my_lev].CellSize();
#if (AMREX_SPACEDIM == 2)
    sum *= dx[0] * dx[1];
#elif (AMREX_SPACEDIM == 3)
    sum *= dx[0] * dx[1] * dx[2];
#endif

    amrex::Print() << "Pred:Sum of mass at time = " << m_cur_time+m_dt << " " << sum << " " << std::endl;
//  amrex::Print() << "Change over time divided by time and area " << (sum - orig_mass) / (m_cur_time+m_dt) / 1.6 <<
//                    " using " << sum << " " << orig_mass << " " << m_cur_time+m_dt << std::endl;
    amrex::Print() << "Pred:Change over time in last time step divided by time and area " << (sum - prev_mass) / m_dt / 1.6 <<
                      " using " << sum << " " << prev_mass << " " << m_dt << std::endl;

#endif

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

#if 1
#ifdef AMREX_USE_EB
    sum = volWgtSum(my_lev,get_density_new_const()[my_lev],0,fact);
#else
    sum = volWgtSum(my_lev,get_density_new_const()[my_lev],0);
#endif

#if (AMREX_SPACEDIM == 2)
    sum *= dx[0] * dx[1];
#elif (AMREX_SPACEDIM == 3)
    sum *= dx[0] * dx[1] * dx[2];
#endif

    amrex::Print() << "Corr:Sum of mass at time = " << m_cur_time+m_dt << " " << sum << " " << std::endl;
//  amrex::Print() << "Change over time divided by time and area " << (sum - orig_mass) / (m_cur_time+m_dt) / 1.6 <<
//                    " using " << sum << " " << orig_mass << " " << m_cur_time+m_dt << std::endl;
    amrex::Print() << "Corr:Change over time in last time step divided by time and area " << (sum - prev_mass) / m_dt / 1.6 <<
                      " using " << sum << " " << prev_mass << " " << m_dt << std::endl;

    prev_mass = sum;
#endif

    // Stop timing current time step
    Real end_step = static_cast<Real>(ParallelDescriptor::second()) - strt_step;
    ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
    if (m_verbose > 0)
    {
        amrex::Print() << "Time per step " << end_step << std::endl;
    }
}

