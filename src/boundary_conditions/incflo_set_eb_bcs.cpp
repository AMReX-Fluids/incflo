#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <incflo.H>

using namespace amrex;


#ifdef AMREX_USE_EB
void
incflo::set_eb_velocity (int lev, amrex::Real time, MultiFab& eb_vel, int nghost)
{
    Geometry const&  gm = Geom(lev);
    const auto& factory = EBFactory(lev, time);

#ifdef AMREX_USE_MOVING_EB
    if ( !eb_vel.ok() ) {
        // We hadn't made the new time factory when leveldata was created, so we
        // define the new time container here
        eb_vel.define(m_leveldata[lev]->velocity.boxArray(), m_leveldata[lev]->velocity.DistributionMap(),
                      AMREX_SPACEDIM, nghost_state(), MFInfo(), factory);
    }
#endif
    eb_vel.setVal(0.);


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(eb_vel, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       const Box& bx = mfi.tilebox();
       const auto& flagfab      = factory.getMultiEBCellFlagFab()[mfi];

       if (flagfab.getType(bx) == FabType::singlevalued) {
          const auto& flags_arr    = flagfab.const_array();
          const auto& eb_vel_arr   = eb_vel[mfi].array();
          const auto& norm_arr     = factory.getBndryNormal()[mfi].const_array();

          bool has_normal = m_eb_flow.has_normal;
          GpuArray<Real, AMREX_SPACEDIM> normal{0.};
          if (has_normal) {
             AMREX_D_TERM(
             normal[0] = m_eb_flow.normal[0];,
             normal[1] = m_eb_flow.normal[1];,
             normal[2] = m_eb_flow.normal[2]);
          }
          Real pad = std::numeric_limits<float>::epsilon();
          Real normal_tol = m_eb_flow.normal_tol;
          Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
          Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

          bool has_comps = false;
          Real eb_vel_mag(0.0);
          GpuArray<amrex::Real,3> eb_vel_comps{0.};
          // Flow is specified as a velocity magnitude
          if ( m_eb_flow.is_mag ) {
              eb_vel_mag = m_eb_flow.vel_mag;
          } else if ( m_eb_flow.is_frequency) {
              has_comps = true;

              const auto& frequency = m_eb_flow.frequency;
              const auto& amplitude = m_eb_flow.amplitude;

              Vector<Real> period(3);
              Vector<Real> amps(3);

              /*
                 --- Adjust frequency to period ---
                 The idea here is that if the frequency is too small, we just neglect it by
                 setting the period to 1 and amplitude to 0.
               */
              if (std::abs(frequency[0]) > 1.e-6){
                  period[0] = 1./frequency[0];
                  amps[0] = amplitude[0];
              } else {
                  period[0] = 1.;
                  amps[0] = 0.;
              }

              if (std::abs(frequency[1]) > 1.e-6){
                  period[1] = 1./frequency[1];
                  amps[1] = amplitude[1];
              } else {
                  period[1] = 1.;
                  amps[1] = 0.;
              }
#if (AMREX_SPACEDIM == 3)
              if (std::abs(frequency[2]) > 1.e-6){
                  period[2] = 1./frequency[2];
                  amps[2] = amplitude[2];
              } else {
                  period[2] = 1.;
                  amps[2] = 0.;
              }
#endif

              Real PI = 3.1415926535897932384;

              AMREX_D_TERM(eb_vel_comps[0] = amps[0] * sin(2*PI*frequency[0]*time) * 2 * PI * frequency[0];,
                           eb_vel_comps[1] = amps[1] * sin(2*PI*frequency[1]*time) * 2 * PI * frequency[1];,
                           eb_vel_comps[2] = amps[2] * sin(2*PI*frequency[2]*time) * 2 * PI * frequency[2];);
          } else {
             has_comps = true;
             const auto& vels = m_eb_flow.velocity;
             AMREX_D_TERM(eb_vel_comps[0] = vels[0];,
                          eb_vel_comps[1] = vels[1];,
                          eb_vel_comps[2] = vels[2]);
          }

          ParallelFor(bx, [flags_arr,eb_vel_arr,norm_arr,has_comps,has_normal,normal,
                 norm_tol_lo, norm_tol_hi,eb_vel_mag,eb_vel_comps]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
             if (flags_arr(i,j,k).isSingleValued()) {
                Real mask = Real(1.0);

                if(has_normal) {
#if (AMREX_SPACEDIM==3)
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1]
                                 + norm_arr(i,j,k,2)*normal[2];
#else
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1];
#endif

                  mask = ((norm_tol_lo <= dotprod) &&
                          (dotprod <= norm_tol_hi)) ? Real(1.0) : Real(0.0);
                }

                if (has_comps) {
                   AMREX_D_TERM(eb_vel_arr(i,j,k,0) = mask*eb_vel_comps[0];,
                                eb_vel_arr(i,j,k,1) = mask*eb_vel_comps[1];,
                                eb_vel_arr(i,j,k,2) = mask*eb_vel_comps[2]);
                } else {
                   // The EB normal points out of the domain so we need to flip the
                   // when using it to convert magnitude to velocity components so
                   // the resulting vector points into the domain.
                   AMREX_D_TERM(eb_vel_arr(i,j,k,0) = -mask*norm_arr(i,j,k,0)*eb_vel_mag;,
                                eb_vel_arr(i,j,k,1) = -mask*norm_arr(i,j,k,1)*eb_vel_mag;,
                                eb_vel_arr(i,j,k,2) = -mask*norm_arr(i,j,k,2)*eb_vel_mag);
                }
             }
           });
       }
     }

     // We make sure to only fill "nghost" ghost cells so we don't accidentally
     // over-write good ghost cell values with unfilled ghost cell values
     IntVect ng_vect(AMREX_D_DECL(nghost,nghost,nghost));
     eb_vel.EnforcePeriodicity(0,AMREX_SPACEDIM,ng_vect,gm.periodicity());
}

void
incflo::set_eb_velocity_for_rotation (int lev, amrex::Real time, MultiFab& eb_vel, int nghost)
{
    Geometry const& gm = Geom(lev);
    const auto dx = gm.CellSizeArray();
    eb_vel.setVal(0.);

    const auto& factory =    EBFactory(lev, time);
    // const auto& factory =
    //    dynamic_cast<EBFArrayBoxFactory const&>(eb_vel.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(eb_vel, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       const Box& bx = mfi.tilebox();
       const auto& flagfab      = factory.getMultiEBCellFlagFab()[mfi];

       if (flagfab.getType(bx) == FabType::singlevalued) {
          const auto& flags_arr    = flagfab.const_array();
          const auto& eb_vel_arr   = eb_vel[mfi].array();
          const auto& bcent_arr     = factory.getBndryCent()[mfi].const_array();

          GpuArray<amrex::Real,3> eb_omega{0.}, eb_cor{0.};
          if ( m_eb_flow.is_omega ) {
             const auto& omega = m_eb_flow.omega;
             eb_omega[0] = omega[0];
             eb_omega[1] = omega[1];
             eb_omega[2] = omega[2];

             const auto& cor = m_eb_flow.center_of_rotation;
             AMREX_D_TERM(eb_cor[0] = cor[0];,
                          eb_cor[1] = cor[1];,
                          eb_cor[2] = cor[2];);
          }

          ParallelFor(bx, [dx,flags_arr,eb_vel_arr,bcent_arr,eb_omega,eb_cor]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
             if (flags_arr(i,j,k).isSingleValued()) {
               Vector<Real> r(3, 0.0);
               RealVect btan;

               // Compute r as vector from EB surface to center of rotation
               AMREX_D_TERM(r[0] = (i + 0.5 + bcent_arr(i,j,k,0))*dx[0] - eb_cor[0];,
                            r[1] = (j + 0.5 + bcent_arr(i,j,k,1))*dx[1] - eb_cor[1];,
                            r[2] = (k + 0.5 + bcent_arr(i,j,k,2))*dx[2] - eb_cor[2];);

               // Compute tangent vector as omega x r
               AMREX_D_TERM(btan[0] = eb_omega[1]*r[2] - eb_omega[2]*r[1];,
                            btan[1] = eb_omega[2]*r[0] - eb_omega[0]*r[2];,
                            btan[2] = eb_omega[0]*r[1] - eb_omega[1]*r[0];);

               AMREX_D_TERM(eb_vel_arr(i,j,k,0) = btan[0];,
                            eb_vel_arr(i,j,k,1) = btan[1];,
                            eb_vel_arr(i,j,k,2) = btan[2];);

             }

           });
       }
     }

     // We make sure to only fill "nghost" ghost cells so we don't accidentally
     // over-write good ghost cell values with unfilled ghost cell values
     IntVect ng_vect(AMREX_D_DECL(nghost,nghost,nghost));
     eb_vel.EnforcePeriodicity(0,AMREX_SPACEDIM,ng_vect,gm.periodicity());
}

void
incflo::set_eb_density (int lev, amrex::Real time, MultiFab& eb_density, int nghost)
{
    Geometry const& gm = Geom(lev);
    eb_density.setVal(0.);

    const auto& factory =    EBFactory(lev, time);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(eb_density, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       const Box& bx = mfi.tilebox();
       const auto& flagfab      = factory.getMultiEBCellFlagFab()[mfi];

       if (flagfab.getType(bx) == FabType::singlevalued) {
          const auto& flags_arr    = flagfab.const_array();
          const auto& eb_density_arr   = eb_density[mfi].array();
          const auto& norm_arr     = factory.getBndryNormal()[mfi].const_array();

          bool has_normal = m_eb_flow.has_normal;
          GpuArray<Real, AMREX_SPACEDIM> normal{0.};
          if (has_normal) {
             AMREX_D_TERM(
             normal[0] = m_eb_flow.normal[0];,
             normal[1] = m_eb_flow.normal[1];,
             normal[2] = m_eb_flow.normal[2]);
          }
          Real pad = std::numeric_limits<float>::epsilon();
          Real normal_tol = m_eb_flow.normal_tol;
          Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
          Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

          Real eb_flow_density = m_eb_flow.density;
// FIXME - this may require that we fill NU cells in the new time arrays also...
          const auto& cc_density = (time == m_cur_time ) ?
              m_leveldata[lev]->density_o[mfi].array() :
              m_leveldata[lev]->density[mfi].array();

          ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
             if (flags_arr(i,j,k).isSingleValued()) {
                Real mask = Real(1.0);

                if(has_normal) {
#if (AMREX_SPACEDIM==3)
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1]
                                 + norm_arr(i,j,k,2)*normal[2];
#else
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1];
#endif

                  mask = ((norm_tol_lo <= dotprod) &&
                          (dotprod <= norm_tol_hi)) ? Real(1.0) : Real(0.0);
                }

#ifdef AMREX_USE_MOVING_EB
                eb_density_arr(i,j,k) = mask*cc_density(i,j,k);
#else
                eb_density_arr(i,j,k) = mask*eb_flow_density;
#endif
             }
           });
       }
     }

     // We make sure to only fill "nghost" ghost cells so we don't accidentally
     // over-write good ghost cell values with unfilled ghost cell values
     IntVect ng_vect(AMREX_D_DECL(nghost,nghost,nghost));
     eb_density.EnforcePeriodicity(0,1,ng_vect,gm.periodicity());
}

void
incflo::set_eb_tracer (int lev, amrex::Real time, MultiFab& eb_tracer, int nghost)
{
    Geometry const& gm = Geom(lev);
    eb_tracer.setVal(0.);

    const auto& factory =    EBFactory(lev, time);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(eb_tracer, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       const Box& bx = mfi.tilebox();
       const auto& flagfab      = factory.getMultiEBCellFlagFab()[mfi];

       if (flagfab.getType(bx) == FabType::singlevalued) {
          const auto& flags_arr     = flagfab.const_array();
          const auto& eb_tracer_arr = eb_tracer[mfi].array();
          const auto& norm_arr      = factory.getBndryNormal()[mfi].const_array();

          bool has_normal = m_eb_flow.has_normal;
          GpuArray<Real, AMREX_SPACEDIM> normal{0.};
          if (has_normal) {
             AMREX_D_TERM(
             normal[0] = m_eb_flow.normal[0];,
             normal[1] = m_eb_flow.normal[1];,
             normal[2] = m_eb_flow.normal[2]);
          }
          Real pad = std::numeric_limits<float>::epsilon();
          Real normal_tol = m_eb_flow.normal_tol;
          Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
          Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

          int num_trac = m_ntrac;
#ifdef AMREX_USE_MOVING_EB
// FIXME - this may require that we fill NU cells in the new time arrays also...
          const auto& cc_tracer = (time == m_cur_time ) ?
              m_leveldata[lev]->tracer_o[mfi].array() :
              m_leveldata[lev]->tracer[mfi].array();
#else
          // Create a device vector
          Gpu::DeviceVector<Real> eb_flow_tracer_dv;
          for(int n(0); n < num_trac; ++n) {
             eb_flow_tracer_dv.push_back(m_eb_flow.tracer[n]);
          }
          Real* eb_flow_tracer = eb_flow_tracer_dv.data();
#endif

          ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
             if (flags_arr(i,j,k).isSingleValued()) {
                Real mask = Real(1.0);

                if(has_normal) {
#if (AMREX_SPACEDIM==3)
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1]
                                 + norm_arr(i,j,k,2)*normal[2];
#else
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1];
#endif

                  mask = ((norm_tol_lo <= dotprod) &&
                          (dotprod <= norm_tol_hi)) ? Real(1.0) : Real(0.0);
                }

                for(int n(0); n<num_trac; n++) {
#ifdef AMREX_USE_MOVING_EB
                    eb_tracer_arr(i,j,k,n) = mask*cc_tracer(i,j,k,n);
#else
                    eb_tracer_arr(i,j,k,n) = mask*eb_flow_tracer[n];
#endif
                }
             }
           });
       }
     }

     // We make sure to only fill "nghost" ghost cells so we don't accidentally
     // over-write good ghost cell values with unfilled ghost cell values
     IntVect ng_vect(AMREX_D_DECL(nghost,nghost,nghost));
     eb_tracer.EnforcePeriodicity(0,m_ntrac,ng_vect,gm.periodicity());
}

#endif
