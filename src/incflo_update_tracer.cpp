#include <incflo.H>

using namespace amrex;

void incflo::update_tracer (StepType step_type, Vector<MultiFab>& tra_forces)
{
    BL_PROFILE("incflo::update_tracer");

    Vector<MultiFab> tra_eta;

    Real new_time = m_cur_time + m_dt;

    if (m_advect_tracer)
    {
        // *************************************************************************************
        // Compute the tracer forcing terms ( forcing for (rho s) if conservative )
        // *************************************************************************************
        for (int lev = 0; lev <= finest_level; ++lev) {
            tra_eta.emplace_back(grids[lev], dmap[lev], m_ntrac, 1, MFInfo(), Factory(lev));
        }

        compute_tracer_diff_coeff(GetVecOfPtrs(tra_eta),1);

        // *************************************************************************************
        // Compute explicit diffusive terms
        // *************************************************************************************
        if (step_type == StepType::Predictor && need_divtau()) {
            compute_laps(get_laps_old(), get_tracer_old_const(), get_density_old_const(),
                         GetVecOfConstPtrs(tra_eta));
        } else if (step_type == StepType::Corrector && m_diff_type == DiffusionType::Explicit) {
            compute_laps(get_laps_new(), get_tracer_new_const(), get_density_new_const(),
                         GetVecOfConstPtrs(tra_eta));
        }

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
