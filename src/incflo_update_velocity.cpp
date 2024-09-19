#include <incflo.H>

using namespace amrex;

void incflo::update_velocity (StepType step_type, Vector<MultiFab>& vel_eta, Vector<MultiFab>& vel_forces)
{
    BL_PROFILE("incflo::update_velocity");

    Real new_time = m_cur_time + m_dt;

    Real l_dt   = m_dt;
    Real l_half = Real(0.5);

    if (step_type == StepType::Predictor) {

        // *************************************************************************************
        // Define (or if advection_type != "MOL", re-define) the forcing terms, without the viscous terms
        //    and using the half-time density
        // *************************************************************************************
        compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_old_const(),
                           get_density_nph_const(), get_tracer_old_const(), get_tracer_new_const());

        // *************************************************************************************
        // Update the velocity
        // *************************************************************************************
        for (int lev = 0; lev <= finest_level; lev++)
        {
        auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& vel = ld.velocity.array(mfi);
            Array4<Real const> const& dvdt = ld.conv_velocity_o.const_array(mfi);
            Array4<Real const> const& vel_f = vel_forces[lev].const_array(mfi);
            Array4<Real const> const& rho_old  = ld.density_o.const_array(mfi);
            Array4<Real const> const& rho_new  = ld.density.const_array(mfi);
            Array4<Real const> const& rho_nph  = ld.density_nph.const_array(mfi);

            if (m_diff_type == DiffusionType::Implicit) {

                if (use_tensor_correction)
                {
                    Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                    // Here divtau_o is the difference of tensor and scalar divtau_o!
                    if (m_advect_momentum) {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                         vel(i,j,k,1) *= rho_old(i,j,k);,
                                         vel(i,j,k,2) *= rho_old(i,j,k););
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                        });
                    }
                } else {
                    if (m_advect_momentum) {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                         vel(i,j,k,1) *= rho_old(i,j,k);,
                                         vel(i,j,k,2) *= rho_old(i,j,k););
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)););
                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)););
                        });
                    }
                }
            }
            else if (m_diff_type == DiffusionType::Crank_Nicolson)
            {

                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                if (m_advect_momentum) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                     vel(i,j,k,1) *= rho_old(i,j,k);,
                                     vel(i,j,k,2) *= rho_old(i,j,k););
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0)+l_half*divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1)+l_half*divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)+l_half*divtau_o(i,j,k,2)););
                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+l_half*divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+l_half*divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+l_half*divtau_o(i,j,k,2)););
                    });
                }
            }
            else if (m_diff_type == DiffusionType::Explicit)
            {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                if (m_advect_momentum) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                     vel(i,j,k,1) *= rho_old(i,j,k);,
                                     vel(i,j,k,2) *= rho_old(i,j,k););
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                    });
                }
            }
        } // mfi
        } // lev

    } else if (step_type == StepType::Corrector) {

        // *************************************************************************************
        // Define the forcing terms to use in the final update (using half-time density)
        // *************************************************************************************
        compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_new_const(),
                           get_density_nph_const(), get_tracer_old_const(), get_tracer_new_const());

        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& vel = ld.velocity.array(mfi);
            Array4<Real const> const& vel_o = ld.velocity_o.const_array(mfi);
            Array4<Real const> const& dvdt = ld.conv_velocity.const_array(mfi);
            Array4<Real const> const& dvdt_o = ld.conv_velocity_o.const_array(mfi);
            Array4<Real const> const& vel_f = vel_forces[lev].const_array(mfi);

            Array4<Real const> const& rho_old  = ld.density_o.const_array(mfi);
            Array4<Real const> const& rho_new  = ld.density.const_array(mfi);
            Array4<Real const> const& rho_nph  = ld.density_nph.const_array(mfi);

            if (m_diff_type == DiffusionType::Explicit)
            {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                Array4<Real const> const& divtau   = ld.divtau.const_array(mfi);

                if (m_advect_momentum) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,0)+  dvdt(i,j,k,0))
                                                  + l_half*(divtau_o(i,j,k,0)+divtau(i,j,k,0))
                                                  + rho_nph(i,j,k) * vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) =  rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                     l_half*(  dvdt_o(i,j,k,1)+  dvdt(i,j,k,1))
                                                   + l_half*(divtau_o(i,j,k,1)+divtau(i,j,k,1))
                                                   + rho_nph(i,j,k) * vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) =  rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                     l_half*(  dvdt_o(i,j,k,2)+  dvdt(i,j,k,2))
                                                   + l_half*(divtau_o(i,j,k,2)+divtau(i,j,k,2))
                                                   + rho_nph(i,j,k) * vel_f(i,j,k,2) ););

                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,0)+  dvdt(i,j,k,0))
                                                  + l_half*(divtau_o(i,j,k,0)+divtau(i,j,k,0))
                                                  + vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,1)+  dvdt(i,j,k,1))
                                                  + l_half*(divtau_o(i,j,k,1)+divtau(i,j,k,1))
                                                  + vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,2)+  dvdt(i,j,k,2))
                                                  + l_half*(divtau_o(i,j,k,2)+divtau(i,j,k,2))
                                                  + vel_f(i,j,k,2) ););
                    });
                }
            }
            else if (m_diff_type == DiffusionType::Crank_Nicolson)
            {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);

                if (m_advect_momentum) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0))
                                                   +l_half*(divtau_o(i,j,k,0) ) + rho_nph(i,j,k)*vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) = rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1))
                                                   +l_half*(divtau_o(i,j,k,1) ) + rho_nph(i,j,k)*vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) = rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2))
                                                   +l_half*(divtau_o(i,j,k,2) ) + rho_nph(i,j,k)*vel_f(i,j,k,2) ););

                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0))
                                                  + l_half* divtau_o(i,j,k,0) + vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1))
                                                  + l_half* divtau_o(i,j,k,1) + vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2))
                                                  + l_half* divtau_o(i,j,k,2) + vel_f(i,j,k,2) ););
                    });
                }
            }
            else if (m_diff_type == DiffusionType::Implicit)
            {
                if (use_tensor_correction)
                {
                    Array4<Real const> const& divtau   = ld.divtau.const_array(mfi);
                    if (m_advect_momentum)
                    {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                        l_half*(dvdt_o(i,j,k,0) +   dvdt(i,j,k,0))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,0) + divtau(i,j,k,0));,
                                         vel(i,j,k,1) = rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                        l_half*(dvdt_o(i,j,k,1) +   dvdt(i,j,k,1))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,1) + divtau(i,j,k,1));,
                                         vel(i,j,k,2) = rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                        l_half*(dvdt_o(i,j,k,2) +   dvdt(i,j,k,2))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,2) + divtau(i,j,k,2)););

                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,0) +   dvdt(i,j,k,0))
                                                      + vel_f(i,j,k,0) + divtau(i,j,k,0));,
                                         vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                        l_half*(dvdt_o(i,j,k,1) +   dvdt(i,j,k,1))
                                                      + vel_f(i,j,k,1) + divtau(i,j,k,1) );,
                                         vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,2) +   dvdt(i,j,k,2))
                                                      + vel_f(i,j,k,2) + divtau(i,j,k,2) ););
                        });
                    }
                } else {
                    if (m_advect_momentum)
                    {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,0) );,
                                         vel(i,j,k,1) = rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,1) );,
                                         vel(i,j,k,2) = rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,2) ););
                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0)) + vel_f(i,j,k,0) );,
                                         vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1)) + vel_f(i,j,k,1) );,
                                         vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2)) + vel_f(i,j,k,2) ););
                        });
                    }
                }
            }
        }
        } // lev
    } // Corrector

    // *************************************************************************************
    // Solve diffusion equation for u* but using eta_old at old time
    // *************************************************************************************
    if (m_diff_type == DiffusionType::Crank_Nicolson || m_diff_type == DiffusionType::Implicit)
    {
        const int ng_diffusion = 1;
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillphysbc_velocity(lev, new_time, m_leveldata[lev]->velocity, ng_diffusion);
            fillphysbc_density (lev, new_time, m_leveldata[lev]->density , ng_diffusion);
        }

        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_dt : l_half*m_dt;
        diffuse_velocity(get_velocity_new(), get_density_new(), GetVecOfConstPtrs(vel_eta), dt_diff);
    }

    // add surface tension
    //fixme: we just consider the surface tension for first tracer
if(0){
    if (m_vof_advect_tracer && m_sigma[0]!=0.){
      VolumeOfFluid*  vof_p = get_volume_of_fluid ();


      for (int lev = 0; lev <= finest_level; lev++)
        {
          auto& ld = *m_leveldata[lev];
          auto const dx = geom[lev].CellSizeArray();

          const auto& ba = ld.density.boxArray();
          const auto& dm = ld.density.DistributionMap();
           const auto& fact = ld.density.Factory();
    Array<MultiFab,AMREX_SPACEDIM> face_val{AMREX_D_DECL(
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                              dm, 1, 0, MFInfo(), fact))};
    MultiFab node_val(amrex::convert(ba,IntVect::TheNodeVector()),dm, AMREX_SPACEDIM, 0 , MFInfo(), fact);

    average_cellcenter_to_face(GetArrOfPtrs(face_val), ld.density, Geom(lev));
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      face_val[idim].invert(m_sigma[0], 0);
      //face_val[idim].setVal(1.0);
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(ld.density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
       // Note nodaltilebox will not include the nodal index beyond boundaries between neighboring
       // titles. Therefore,if we want to use face values (i.e., face_val) immediately below (commented
       // out), we must create index space for the face-centered values of the tiled region
       // (i.e., surroundingNodes()).Since we currently calculate all face values in all boxes and then
       // convert them to node-centered value, it is good that we use (nodaltilebox()) to avoid the
       // repeated calculation of face-centered values at cell faces which are shared by two tiles.
       Box const& xbx = mfi.nodaltilebox(0);
       Box const& ybx = mfi.nodaltilebox(1);
       Box const& zbx = mfi.nodaltilebox(2);
       Array4<Real const> const& tra   = ld.tracer.const_array(mfi);
       Array4<Real const> const& kap   = vof_p->kappa[lev].const_array(mfi);
       AMREX_D_TERM( Array4<Real > const& xfv = face_val[0].array(mfi);,
                     Array4<Real > const& yfv = face_val[1].array(mfi);,
                     Array4<Real > const& zfv = face_val[2].array(mfi););

       ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         Real kaf;
         if(kap(i,j,k,0)!=VOF_NODATA && kap(i-1,j,k,0)!=VOF_NODATA)
            kaf=Real(0.5)*(kap(i,j,k,0)+kap(i-1,j,k,0));
         else if (kap(i,j,k,0)!=VOF_NODATA)
            kaf=kap(i,j,k);
         else if (kap(i-1,j,k,0)!=VOF_NODATA)
            kaf=kap(i-1,j,k);
         else
            kaf=0.;
          xfv(i,j,k) *= kaf*(tra(i,j,k,0)-tra(i-1,j,k,0))/dx[0];
       });

       ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         Real kaf;
         if(kap(i,j,k,0)!=VOF_NODATA && kap(i,j-1,k,0)!=VOF_NODATA)
            kaf=Real(0.5)*(kap(i,j,k,0)+kap(i,j-1,k,0));
         else if (kap(i,j,k,0)!=VOF_NODATA)
            kaf=kap(i,j,k);
         else if (kap(i,j-1,k,0)!=VOF_NODATA)
            kaf=kap(i,j-1,k);
         else
            kaf=0.;
          yfv(i,j,k) *= kaf*(tra(i,j,k,0)-tra(i,j-1,k,0))/dx[1];
       });
#if AMREX_SPACEDIM == 3
       ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         Real kaf;
         if(kap(i,j,k,0)!=VOF_NODATA && kap(i,j,k-1,0)!=VOF_NODATA)
            kaf=Real(0.5)*(kap(i,j,k,0)+kap(i,j,k-1,0));
         else if (kap(i,j,k,0)!=VOF_NODATA)
            kaf=kap(i,j,k);
         else if (kap(i,j,k-1,0)!=VOF_NODATA)
            kaf=kap(i,j,k-1);
         else
            kaf=0.;
          zfv(i,j,k) *= kaf*(tra(i,j,k,0)-tra(i,j,k-1,0))/dx[2];

     /*    if(i==8&&j==2&&k==8){
      Print()<<"zbx   "<<"low  "<<tra(i,j,k,0)<<" high " <<tra(i,j,k-1,0)<<"  "
             << kaf<<"  "<<zfv(i,j,k)<<"\n";

         }    */
       });
#endif
    }
static int oct[3][2] = { { 1, 2 }, { 0, 2 }, { 0, 1 } };
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(vel_forces[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
       Box const& bx  = mfi.tilebox();
       //We will immediately use the nodal values so we must use surroundingNodes()instead of nodaltilebox()
       //See the previous comments for calculating face-centered value.
       Box const& nbx = surroundingNodes(bx);
       //Box const& nbx = mfi.nodaltilebox();
       Box const& vbx  = mfi.validbox();
       Box const& xvbx = surroundingNodes(vbx,0);
       Box const& yvbx = surroundingNodes(vbx,1);
       Box const& zvbx = surroundingNodes(vbx,2);
       AMREX_D_TERM( Array4<Real > const& xfv = face_val[0].array(mfi);,
                     Array4<Real > const& yfv = face_val[1].array(mfi);,
                     Array4<Real > const& zfv = face_val[2].array(mfi););
       Array4<Real > const& nv = node_val.array(mfi);
       ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         //nv(i,j,k,0)=0.25*(xfv(i,j,k,0)+xfv(i,j-1,k,0)+xfv(i,j,k-1,0)+xfv(i,j-1,k-1,0));
         //nv(i,j,k,1)=0.25*(yfv(i,j,k,0)+yfv(i,j,k-1,0)+yfv(i-1,j,k,0)+yfv(i-1,j,k-1,0));
         //nv(i,j,k,2)=0.25*(zfv(i,j,k,0)+zfv(i-1,j,k,0)+zfv(i,j-1,k,0)+zfv(i-1,j-1,k,0));

         for (int dim = 0; dim < AMREX_SPACEDIM; ++dim){
           nv(i,j,k,dim)=0.;
           int nt=0;
           for (int nn=0; nn<4;++nn){
             Array<int,3> in{i, j, k};
             if(nn==1)
               in[oct[dim][0]]-=1;
             else if (nn==2)
               in[oct[dim][1]]-=1;
             else if (nn==3) {
               in[oct[dim][0]]-=1;
               in[oct[dim][1]]-=1;
             }
             if (dim==0&& xvbx.contains(in[0],in[1],in[2])){
                nv(i,j,k,dim)+= xfv(in[0],in[1],in[2]);
                nt++;
             }
             else if (dim==1&& yvbx.contains(in[0],in[1],in[2])){
                nv(i,j,k,dim)+= yfv(in[0],in[1],in[2]);
                nt++;
             }
#if AMREX_SPACEDIM == 3
             else if (dim==2&& zvbx.contains(in[0],in[1],in[2])){
                nv(i,j,k,dim)+= zfv(in[0],in[1],in[2]);
                nt++;
             }
#endif
           }
           nv(i,j,k,dim)/= nt;

         }
       });

       Array4<Real>  const& forarr  = vof_p->force[lev].array(mfi);
       Array4<Real> const& vel = ld.velocity.array(mfi);
       ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         for (int dim = 0; dim < AMREX_SPACEDIM; ++dim){
           vel(i,j,k,dim) -= Real(0.125)*(nv(i,j  ,k  ,dim) + nv(i+1,j  ,k  ,dim)
                                          + nv(i,j+1,k  ,dim) + nv(i+1,j+1,k  ,dim)
                                          + nv(i,j  ,k+1,dim) + nv(i+1,j  ,k+1,dim)
                                          + nv(i,j+1,k+1,dim) + nv(i+1,j+1,k+1,dim))*m_dt;
           forarr(i,j,k,dim) = -Real(0.125)*(nv(i,j  ,k  ,dim) + nv(i+1,j  ,k  ,dim)
                                          + nv(i,j+1,k  ,dim) + nv(i+1,j+1,k  ,dim)
                                          + nv(i,j  ,k+1,dim) + nv(i+1,j  ,k+1,dim)
                                          + nv(i,j+1,k+1,dim) + nv(i+1,j+1,k+1,dim));
         }

       });


    }




       } // end lev

    }
 }


}

