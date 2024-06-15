#include <incflo.H>

using namespace amrex;

void incflo::update_density (StepType step_type)
{
    BL_PROFILE("incflo::update_density");

    int ng = (step_type == StepType::Corrector) ? 0 : 1;

    Real l_dt = m_dt;

    if (m_constant_density)
    {
        for (int lev = 0; lev <= finest_level; lev++)
            MultiFab::Copy(m_leveldata[lev]->density_nph, m_leveldata[lev]->density_o, 0, 0, 1, ng);

    } else {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real const> const& rho_old  = ld.density_o.const_array(mfi);
                Array4<Real      > const& rho_new  = ld.density.array(mfi);
                Array4<Real      > const& rho_nph  = ld.density_nph.array(mfi);
                Array4<Real const> const& drdt_o = ld.conv_density_o.const_array(mfi);

                if (step_type == StepType::Predictor) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        rho_new(i,j,k) = rho_old(i,j,k) + l_dt * drdt_o(i,j,k);
                    });

                } else if (step_type == StepType::Corrector) {
                    Array4<Real const> const& drdt_n = ld.conv_density.const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        rho_new(i,j,k) = rho_old(i,j,k) + l_dt * Real(0.5) * (drdt_o(i,j,k) + drdt_n(i,j,k));
                        rho_nph(i,j,k) = Real(0.5) * (rho_old(i,j,k) + rho_new(i,j,k));
                    });
                }
            } // mfi

            if (step_type == StepType::Predictor) {
                // Fill ghost cells of the new density field so that we can define density_nph
                //      on the valid region grown by 1
                fillpatch_density(lev, m_t_new[lev], ld.density, ng);

                for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    Box const& gbx = mfi.growntilebox(1);
                    Array4<Real  const> const& rho_new  = ld.density.const_array(mfi);
                    Array4<Real  const> const& rho_old  = ld.density_o.const_array(mfi);
                    Array4<Real       > const& rho_nph  = ld.density_nph.array(mfi);

                    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        rho_nph(i,j,k) = Real(0.5) * (rho_old(i,j,k) + rho_new(i,j,k));
                    });
                } // mfi
            } // Predictor
        } // lev

        // Average down solution
        for (int lev = finest_level-1; lev >= 0; --lev) {
#ifdef AMREX_USE_EB
            amrex::EB_average_down(m_leveldata[lev+1]->density, m_leveldata[lev]->density,
                                   0, 1, refRatio(lev));
#else
            amrex::average_down(m_leveldata[lev+1]->density, m_leveldata[lev]->density,
                                0, 1, refRatio(lev));
#endif
        }
    } // not constant density
}
