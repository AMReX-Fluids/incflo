#include <incflo.H>

using namespace amrex;

void incflo::update_temperature (StepType step_type, Vector<MultiFab>& tem_eta, Vector<MultiFab>& scratch)
{
    BL_PROFILE("incflo::update_temperature");

    if (!m_use_temperature) { return; }

    Real const  new_time = m_cur_time + m_dt;
    Real const half_time = m_cur_time + m_dt/2.;

    // *************************************************************************************
    // Compute the temperature forcing terms
    // *************************************************************************************
    compute_tem_forces(half_time, GetVecOfPtrs(scratch));

    // *************************************************************************************
    // Compute explicit diffusive term (if corrector)
    // *************************************************************************************
    if (step_type == StepType::Corrector)
    {
        compute_temperature_diff_coeff(new_time, GetVecOfPtrs(tem_eta));
        if (m_diff_type == DiffusionType::Explicit) {
            compute_laps_T(get_laps_new(), get_temperature_new_const(), GetVecOfConstPtrs(tem_eta));
        }
    }

    // *************************************************************************************
    // Update the temperature with time-explicit terms
    // *************************************************************************************
    if (step_type == StepType::Predictor) {
        constexpr Real m_half = Real(0.5);
        Real l_dt = m_dt;

        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real const> const& tem_o   = ld.temperature_o.const_array(mfi);
                Array4<Real      > const& tem     = ld.temperature.array(mfi);
                Array4<Real const> const& rho_h   = ld.density_nph.const_array(mfi);
                Array4<Real const> const& dtdt_o  = ld.conv_temperature_o.const_array(mfi);
                // temperature forcing term (Q) is in scratch
                Array4<Real      > const& tem_f   = scratch[lev].array(mfi);

                FArrayBox cp_fab(bx, 1, The_Async_Arena());
                compute_cp(lev, mfi, cp_fab);
                Array4<Real      > const& cp      = cp_fab.array();

                if (m_diff_type == DiffusionType::Explicit)
                {
                    Array4<Real const> const& laps_o = ld.laps_tem_o.const_array(mfi);

                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        tem(i,j,k) = tem_o(i,j,k) + l_dt *
                            ( dtdt_o(i,j,k) + (tem_f(i,j,k) + laps_o(i,j,k))/(rho_h(i,j,k) * cp(i,j,k)) );
                    });
                }
                else if (m_diff_type == DiffusionType::Crank_Nicolson)
                {
                    Array4<Real const> const& laps_o = ld.laps_tem_o.const_array(mfi);

                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        tem(i,j,k) = tem_o(i,j,k) + l_dt *
                            ( dtdt_o(i,j,k) + (tem_f(i,j,k) + m_half*laps_o(i,j,k))/(rho_h(i,j,k) * cp(i,j,k)) );
                        // Save rhoCp for use in implicit solve.
                        // Reuse scratch space since we are done with forcing now.
                        tem_f(i,j,k) = rho_h(i,j,k) * cp(i,j,k);
                    });
                }
                else if (m_diff_type == DiffusionType::Implicit)
                {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        tem(i,j,k) = tem_o(i,j,k) + l_dt *
                            (dtdt_o(i,j,k) + tem_f(i,j,k)) / (rho_h(i,j,k) * cp(i,j,k));
                        // Save rhoCp for use in implicit solve.
                        // Reuse scratch space since we are done with forcing now.
                        tem_f(i,j,k) = rho_h(i,j,k) * cp(i,j,k);
                    });
                }
            } // mfi
        } // lev

    } else if (step_type == StepType::Corrector) {
        Abort("incflo::update_temperature does not yet work with the corrector");
    }

    // *************************************************************************************
    // Solve implicit diffusion equation for temperature
    // *************************************************************************************
    if (m_diff_type == DiffusionType::Crank_Nicolson || m_diff_type == DiffusionType::Implicit)
    {
        const int ng_diffusion = 1;
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillphysbc_temperature(lev, new_time, m_leveldata[lev]->temperature, ng_diffusion);
        }
        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_dt : Real(0.5)*m_dt;
        // scratch holds rhoCp
        diffuse_temperature(get_temperature_new(), GetVecOfPtrs(scratch), GetVecOfConstPtrs(tem_eta),
                            dt_diff);
    }
    else
    {
        // Need to average down temperature since the diffusion solver didn't do it for us.
        for (int lev = finest_level-1; lev >= 0; --lev) {
#ifdef AMREX_USE_EB
            amrex::EB_average_down(m_leveldata[lev+1]->temperature, m_leveldata[lev]->temperature,
                                   0, m_ntrac, refRatio(lev));
#else
            amrex::average_down(m_leveldata[lev+1]->temperature, m_leveldata[lev]->temperature,
                                0, m_ntrac, refRatio(lev));
#endif
        }
    }
}
