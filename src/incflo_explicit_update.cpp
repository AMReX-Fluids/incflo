#include <incflo.H>

using namespace amrex;

void incflo::tracer_explicit_update (Vector<MultiFab> const& tra_forces)
{
    if (m_advect_tracer == 0) { return; }

    constexpr Real m_half = Real(0.5);
    Real l_dt = m_dt;
    int l_ntrac = m_ntrac;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real const> const& tra_o   = ld.tracer_o.const_array(mfi);
            Array4<Real const> const& rho_o   = ld.density_o.const_array(mfi);
            Array4<Real> const& tra           = ld.tracer.array(mfi);
            Array4<Real const> const& rho     = ld.density.const_array(mfi);
            Array4<Real const> const& dtdt_o  = ld.conv_tracer_o.const_array(mfi);
            Array4<Real const> const& tra_f   = tra_forces[lev].const_array(mfi);

            auto const* iconserv = get_tracer_iconserv_device_ptr();

            if (m_diff_type == DiffusionType::Explicit)
            {
                Array4<Real const> const& laps_o = ld.laps_o.const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // If conservative
                    // (rho trac)^new = (rho trac)^old + dt * (
                    //                   div(rho trac u) + div (mu grad trac) + rho * f_t )
                    // else non-conservative
                    // (trac)^new = (trac)^old + dt * (
                    //               u dot grad trac + div (mu grad trac) + f_t )
                    for (int n = 0; n < l_ntrac; ++n)
                    {
                        if ( iconserv[n] ) {
                            Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) + laps_o(i,j,k,n) );
                            tra(i,j,k,n) = tra_new/rho(i,j,k);
                        } else {
                            tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) + laps_o(i,j,k,n) );
                        }
                    }
                });
            }
            else if (m_diff_type == DiffusionType::Crank_Nicolson)
            {
                Array4<Real const> const& laps_o = ld.laps_o.const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int n = 0; n < l_ntrac; ++n)
                    {
                        if ( iconserv[n] ) {
                            Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) + m_half * laps_o(i,j,k,n) );
                            tra(i,j,k,n) = tra_new/rho(i,j,k);
                        } else {
                            tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) + m_half * laps_o(i,j,k,n) );
                        }
                    }
                });
            }
            else if (m_diff_type == DiffusionType::Implicit)
            {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int n = 0; n < l_ntrac; ++n)
                    {
                        if ( iconserv[n] ) {
                            Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) );
                            tra(i,j,k,n) = tra_new/rho(i,j,k);
                        } else {
                            tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) );
                        }

                    }
                });
            }
        } // mfi
    } // lev
}

void incflo::tracer_explicit_update_corrector (Vector<MultiFab> const& tra_forces)
{
    if (m_advect_tracer == 0) { return; }

    constexpr Real m_half = Real(0.5);
    Real l_dt = m_dt;
    int l_ntrac = m_ntrac;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real const> const& tra_o   = ld.tracer_o.const_array(mfi);
            Array4<Real const> const& rho_o   = ld.density_o.const_array(mfi);
            Array4<Real      > const& tra     = ld.tracer.array(mfi);
            Array4<Real const> const& rho     = ld.density.const_array(mfi);
            Array4<Real const> const& dtdt_o  = ld.conv_tracer_o.const_array(mfi);
            Array4<Real const> const& dtdt    = ld.conv_tracer.const_array(mfi);
            Array4<Real const> const& tra_f   = tra_forces[lev].const_array(mfi);
            auto const* iconserv = get_tracer_iconserv_device_ptr();

            if (m_diff_type == DiffusionType::Explicit)
            {
                Array4<Real const> const& laps_o = ld.laps_o.const_array(mfi);
                Array4<Real const> const& laps   = ld.laps.const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int n = 0; n < l_ntrac; ++n)
                    {
                        if ( iconserv[n] ) {
                            Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt * (
                                m_half*(  dtdt(i,j,k,n) + dtdt_o(i,j,k,n))
                                +m_half*(laps_o(i,j,k,n) +   laps(i,j,k,n))
                                +    tra_f(i,j,k,n) );

                            tra(i,j,k,n) = tra_new / rho(i,j,k);
                        } else {
                            tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt * (
                                m_half*(  dtdt(i,j,k,n) + dtdt_o(i,j,k,n))
                                +m_half*(laps_o(i,j,k,n) +   laps(i,j,k,n))
                                +    tra_f(i,j,k,n) );
                        }
                    }
                });
            }
            else if (m_diff_type == DiffusionType::Crank_Nicolson)
            {
                Array4<Real const> const& laps_o = ld.laps_o.const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int n = 0; n < l_ntrac; ++n)
                    {
                        if ( iconserv[n] ) {
                            Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt * (
                                m_half*(  dtdt(i,j,k,n) + dtdt_o(i,j,k,n))
                                +m_half*(laps_o(i,j,k,n)                  )
                                +    tra_f(i,j,k,n) );

                            tra(i,j,k,n) = tra_new / rho(i,j,k);
                        } else {
                            tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt * (
                                m_half*(  dtdt(i,j,k,n) + dtdt_o(i,j,k,n))
                                +m_half*(laps_o(i,j,k,n)                  )
                                +    tra_f(i,j,k,n) );
                        }
                    }
                });
            }
            else if (m_diff_type == DiffusionType::Implicit)
            {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int n = 0; n < l_ntrac; ++n)
                    {
                        if ( iconserv[n] ) {
                            Real tra_new = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt * (
                                m_half*(  dtdt(i,j,k,n)+dtdt_o(i,j,k,n))
                                +      tra_f(i,j,k,n) );

                            tra(i,j,k,n) = tra_new / rho(i,j,k);
                        } else {
                            tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt * (
                                m_half*(  dtdt(i,j,k,n)+dtdt_o(i,j,k,n))
                                +      tra_f(i,j,k,n) );
                        }
                    }
                });
            }
        } // mfi
    } // lev
}
