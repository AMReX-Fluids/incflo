#include <incflo.H>

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

    copy_from_new_to_old_velocity();
    copy_from_new_to_old_density();
    copy_from_new_to_old_tracer();

    
    int ng = nghost_state();
    for (int lev = 0; lev <= finest_level; ++lev) {
        fillpatch_velocity(lev, m_t_old[lev], m_leveldata[lev]->velocity_o, ng);
        fillpatch_density(lev, m_t_old[lev], m_leveldata[lev]->density_o, ng);
        if (m_advect_tracer) {
            fillpatch_tracer(lev, m_t_old[lev], m_leveldata[lev]->tracer_o, ng);
        }
    }

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
    // **********************************************************************************************
    //
    // Update the moving geometry and arrays
    //
    // **********************************************************************************************

    VisMF::Write(m_leveldata[0]->density_o,"do1");
//    if (!incremental_projection) {
    for (int lev = 0; lev <= finest_level; lev++)
    {
	MakeNewGeometry(lev,m_t_new[lev]);
//    }
// Need to fill the to be NU cells here
     // Or we stick with special treatment of NU cells in apply and redistribute.cpp
// Need to update the EBFactory for this to work...
    
	// Now let's make sure to fill cells that were previously covered but are now cut cell
	// Not sure we need to do the new MFs, maybe could get by with just the olds
	//amrex::Print() << "Fill Velocity" << std::endl;
	EB_fill_uncovered(lev, m_leveldata[lev]->velocity  , m_leveldata[lev]->velocity  );
	EB_fill_uncovered(lev, m_leveldata[lev]->velocity_o, m_leveldata[lev]->velocity_o);

	//amrex::Print() << "\nFill density" << std::endl;
	EB_fill_uncovered(lev, m_leveldata[lev]->density   , m_leveldata[lev]->density   );
	EB_fill_uncovered(lev, m_leveldata[lev]->density_o , m_leveldata[lev]->density_o );

	if (m_ntrac > 0) {
	    EB_fill_uncovered(lev, m_leveldata[lev]->tracer   , m_leveldata[lev]->tracer  );
	    EB_fill_uncovered(lev, m_leveldata[lev]->tracer_o , m_leveldata[lev]->tracer_o);
	}

	//FIXME - will need to be more careful here when adding diffusion...
	// problem here in that divtau hasn't been computed yet...
	// EB_fill_uncovered_with_zero(lev, m_leveldata[lev]->divtau_o, m_leveldata[lev]->divtau_o);
	Print()<<"Setting divtau to zero..."<<std::endl;       
	m_leveldata[lev]->divtau_o.setVal(0.0);

// FIXME also need to worry about all the pieces of U*, forces, etc...

	VisMF::Write(m_leveldata[lev]->gp,"gp");
	//amrex::Print() << "\nFill gp" << std::endl;
	EB_fill_uncovered(lev, m_leveldata[lev]->gp      , m_leveldata[lev]->gp  );

	// This function is for cell-centered data. Not garaunteed to be correct for
	// nodal or face centered data...
	//amrex::Print() << "\nFill p_nd" << std::endl;
	EB_fill_uncovered(lev, m_leveldata[lev]->p_nd    , m_leveldata[lev]->p_nd);
    }
#endif
    VisMF::Write(m_leveldata[0]->velocity,"vel");
    VisMF::Write(m_leveldata[0]->velocity_o,"velo");
    VisMF::Write(m_leveldata[0]->density_o,"do2");
    
    ApplyPredictor();

//FIXME
    // this will overwrite the previous time plotfile
    //WritePlotFile();
    static int count=0; count++;
    if (count>2) Abort();
    
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

