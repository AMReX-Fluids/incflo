#include <incflo.H>

using namespace amrex;

void incflo::update_tracer (StepType step_type, Vector<MultiFab>& tra_eta, Vector<MultiFab>& tra_forces)
{
    BL_PROFILE("incflo::update_tracer");

    Real new_time = m_cur_time + m_dt;

    if (m_advect_tracer)
    {
        // *************************************************************************************
        // Compute diffusion coefficient for tracer
        // *************************************************************************************

        // *************************************************************************************
        // Compute the tracer forcing terms (forcing for (rho s), not for s)
        // *************************************************************************************
        compute_tra_forces(GetVecOfPtrs(tra_forces),  get_density_nph_const());

        // *************************************************************************************
        // Compute explicit diffusive term (if corrector)
        // *************************************************************************************
        if (step_type == StepType::Corrector)
        {
            compute_tracer_diff_coeff(GetVecOfPtrs(tra_eta),1);
            if (m_diff_type == DiffusionType::Explicit) {
                compute_laps(get_laps_new(), get_tracer_new_const(), get_density_new_const(),
                             GetVecOfConstPtrs(tra_eta));
            }
        }

        // *************************************************************************************
        // Update the tracer next (note that dtdt already has rho in it)
        // (rho trac)^new = (rho trac)^old + dt * (
        //                   div(rho trac u) + div (mu grad trac) + rho * f_t
        // *************************************************************************************
        if (step_type == StepType::Predictor) {
            tracer_explicit_update(tra_forces);
        } else if (step_type == StepType::Corrector) {
            tracer_explicit_update_corrector(tra_forces);
        }

        // *************************************************************************************
        // Solve diffusion equation for tracer
        // *************************************************************************************
        if (m_diff_type == DiffusionType::Crank_Nicolson || m_diff_type == DiffusionType::Implicit)
        {
            const int ng_diffusion = 1;
            for (int lev = 0; lev <= finest_level; ++lev)
                fillphysbc_tracer(lev, new_time, m_leveldata[lev]->tracer, ng_diffusion);

            Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_dt : Real(0.5)*m_dt;
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
    } // advect tracer
}
