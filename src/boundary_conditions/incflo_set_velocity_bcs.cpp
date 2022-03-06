#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <incflo.H>

using namespace amrex;

void
incflo::set_inflow_velocity (int lev, amrex::Real time, MultiFab& vel, int nghost)
{
    Geometry const& gm = Geom(lev);
    Box const& domain = gm.growPeriodicDomain(nghost);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        Orientation olo(dir,Orientation::low);
        Orientation ohi(dir,Orientation::high);
        if (m_bc_type[olo] == BC::mass_inflow || m_bc_type[ohi] == BC::mass_inflow) {
            Box dlo = (m_bc_type[olo] == BC::mass_inflow) ? amrex::adjCellLo(domain,dir,nghost) : Box();
            Box dhi = (m_bc_type[ohi] == BC::mass_inflow) ? amrex::adjCellHi(domain,dir,nghost) : Box();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
                Box const& gbx = amrex::grow(mfi.validbox(),nghost);
                Box blo = gbx & dlo;
                Box bhi = gbx & dhi;
                Array4<Real> const& v = vel[mfi].array();
                int gid = mfi.index();
                if (blo.ok()) {
                    prob_set_inflow_velocity(gid, olo, blo, v, lev, time);
                }
                if (bhi.ok()) {
                    prob_set_inflow_velocity(gid, ohi, bhi, v, lev, time);
                }
            }
        }
    }
    // We make sure to only fill "nghost" ghost cells so we don't accidentally 
    // over-write good ghost cell values with unfilled ghost cell values 
    IntVect ng_vect(AMREX_D_DECL(nghost,nghost,nghost));
    vel.EnforcePeriodicity(0,AMREX_SPACEDIM,ng_vect,gm.periodicity());
}

void
incflo::set_eb_velocity (int lev, amrex::Real time, MultiFab& eb_vel, int nghost, 
      const amrex::Real eb_vel_mag)
{
    Geometry const& gm = Geom(lev);
    eb_vel.setVal(0.);

    const auto& factory =
       dynamic_cast<EBFArrayBoxFactory const&>(eb_vel.Factory());

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

          AMREX_D_TERM(
                const auto& apx = factory.getAreaFrac()[0]->const_array(mfi);,
                const auto& apy = factory.getAreaFrac()[1]->const_array(mfi);,
                const auto& apz = factory.getAreaFrac()[2]->const_array(mfi));

           ParallelFor(bx, [flags_arr,eb_vel_arr,norm_arr,
                 eb_vel_mag,AMREX_D_DECL(apx,apy,apz)]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
             if (flags_arr(i,j,k).isSingleValued()) {
               AMREX_D_TERM(
                  Real apxm = apx(i  ,j  ,k  );,
                  Real apym = apy(i  ,j  ,k  );,
                  Real apzm = apz(i  ,j  ,k  ));

               AMREX_D_TERM(
                  Real apxp = apx(i+1,j  ,k  );,
                  Real apyp = apy(i  ,j+1,k  );,
                  Real apzp = apz(i  ,j  ,k+1));

               AMREX_D_TERM(
                  Real dapx = apxm-apxp;,
                  Real dapy = apym-apyp;,
                  Real dapz = apzm-apzp);

#if (AMREX_SPACEDIM == 3)
               Real anorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
#else
               Real anorm = std::sqrt(dapx*dapx+dapy*dapy);
#endif
               Real anorminv = 1.0/anorm;

               AMREX_D_TERM(
                     Real anrmx = dapx * anorminv;,
                     Real anrmy = dapy * anorminv;,
                     Real anrmz = dapz * anorminv);

               AMREX_D_TERM(
                     eb_vel_arr(i,j,k,0) = -anrmx*eb_vel_mag;,
                     eb_vel_arr(i,j,k,1) = -anrmy*eb_vel_mag;,
                     eb_vel_arr(i,j,k,2) = -anrmz*eb_vel_mag);
             }
           });
        }
     }

    // We make sure to only fill "nghost" ghost cells so we don't accidentally 
    // over-write good ghost cell values with unfilled ghost cell values 
    IntVect ng_vect(AMREX_D_DECL(nghost,nghost,nghost));
    eb_vel.EnforcePeriodicity(0,AMREX_SPACEDIM,ng_vect,gm.periodicity());
}
