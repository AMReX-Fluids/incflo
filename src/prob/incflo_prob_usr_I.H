#ifndef INCFLO_PROB_USR_I_H_
#define INCFLO_PROB_USR_I_H_

// This header is include by incflo.H.  These functions are members of incflo class.

    void init_heated_ground (amrex::Box const& vbx, amrex::Box const& gbx,
                             amrex::Array4<amrex::Real> const& p,
                             amrex::Array4<amrex::Real> const& vel,
                             amrex::Array4<amrex::Real> const& density,
                             amrex::Array4<amrex::Real> const& tracer,
                             amrex::Box const& domain,
                             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_taylor_green (amrex::Box const& vbx, amrex::Box const& gbx,
                            amrex::Array4<amrex::Real> const& p,
                            amrex::Array4<amrex::Real> const& vel,
                            amrex::Array4<amrex::Real> const& density,
                            amrex::Array4<amrex::Real> const& tracer,
                            amrex::Box const& domain,
                            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_taylor_green3d (amrex::Box const& vbx, amrex::Box const& gbx,
                              amrex::Array4<amrex::Real> const& p,
                              amrex::Array4<amrex::Real> const& vel,
                              amrex::Array4<amrex::Real> const& density,
                              amrex::Array4<amrex::Real> const& tracer,
                              amrex::Box const& domain,
                              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_taylor_vortex (amrex::Box const& vbx, amrex::Box const& gbx,
                             amrex::Array4<amrex::Real> const& p,
                             amrex::Array4<amrex::Real> const& vel,
                             amrex::Array4<amrex::Real> const& density,
                             amrex::Array4<amrex::Real> const& tracer,
                             amrex::Box const& domain,
                             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_couette (amrex::Box const& vbx, amrex::Box const& gbx,
                       amrex::Array4<amrex::Real> const& p,
                       amrex::Array4<amrex::Real> const& vel,
                       amrex::Array4<amrex::Real> const& density,
                       amrex::Array4<amrex::Real> const& tracer,
                       amrex::Box const& domain,
                       amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                       amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                       amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_channel_slant (amrex::Box const& vbx, amrex::Box const& gbx,
                             amrex::Array4<amrex::Real> const& p,
                             amrex::Array4<amrex::Real> const& vel,
                             amrex::Array4<amrex::Real> const& density,
                             amrex::Array4<amrex::Real> const& tracer,
                             amrex::Box const& domain,
                             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_rayleigh_taylor (amrex::Box const& vbx, amrex::Box const& gbx,
                               amrex::Array4<amrex::Real> const& p,
                               amrex::Array4<amrex::Real> const& vel,
                               amrex::Array4<amrex::Real> const& density,
                               amrex::Array4<amrex::Real> const& tracer,
                               amrex::Box const& domain,
                               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_tuscan (amrex::Box const& vbx, amrex::Box const& gbx,
                      amrex::Array4<amrex::Real> const& p,
                      amrex::Array4<amrex::Real> const& vel,
                      amrex::Array4<amrex::Real> const& density,
                      amrex::Array4<amrex::Real> const& tracer,
                      amrex::Box const& domain,
                      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_boussinesq_bubble (amrex::Box const& vbx, amrex::Box const& gbx,
                                 amrex::Array4<amrex::Real> const& p,
                                 amrex::Array4<amrex::Real> const& vel,
                                 amrex::Array4<amrex::Real> const& density,
                                 amrex::Array4<amrex::Real> const& tracer,
                                 amrex::Box const& domain,
                                 amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                                 amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                                 amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_periodic_tracer (amrex::Box const& vbx, amrex::Box const& gbx,
                               amrex::Array4<amrex::Real> const& p,
                               amrex::Array4<amrex::Real> const& vel,
                               amrex::Array4<amrex::Real> const& density,
                               amrex::Array4<amrex::Real> const& tracer,
                               amrex::Box const& domain,
                               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_double_shear_layer (amrex::Box const& vbx, amrex::Box const& gbx,
                                  amrex::Array4<amrex::Real> const& p,
                                  amrex::Array4<amrex::Real> const& vel,
                                  amrex::Array4<amrex::Real> const& density,
                                  amrex::Array4<amrex::Real> const& tracer,
                                  amrex::Box const& domain,
                                  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                                  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                                  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_plane_poiseuille (amrex::Box const& vbx, amrex::Box const& gbx,
                                amrex::Array4<amrex::Real> const& p,
                                amrex::Array4<amrex::Real> const& vel,
                                amrex::Array4<amrex::Real> const& density,
                                amrex::Array4<amrex::Real> const& tracer,
                                amrex::Box const& domain,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);
    void init_flow_in_box (amrex::Box const& vbx, amrex::Box const& gbx,
                           amrex::Array4<amrex::Real> const& p,
                           amrex::Array4<amrex::Real> const& vel,
                           amrex::Array4<amrex::Real> const& density,
                           amrex::Array4<amrex::Real> const& tracer,
                           amrex::Box const& domain,
                           amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                           amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                           amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);
#endif
