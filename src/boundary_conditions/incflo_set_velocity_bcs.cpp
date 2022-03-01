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
incflo::set_eb_velocity (int lev, amrex::Real time, MultiFab& vel, int nghost, 
      const Vector<Real>& eb_velocity, const Vector<Real>& eb_normal)
{
    Geometry const& gm = Geom(lev);
    vel.setVal(0.);

    const auto& factory =
       dynamic_cast<EBFArrayBoxFactory const&>(vel.Factory());


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(vel, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       const Box& bx = mfi.tilebox();

       const auto& vel_arr      = vel[mfi].array();
       const auto& flags_arr    = factory.getMultiEBCellFlagFab()[mfi].const_array();

       const auto& norm_arr  = factory.getBndryNormal()[mfi].const_array();

        ParallelFor(bx, [flags_arr,vel_arr,norm_arr,eb_velocity,eb_normal]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (flags_arr(i,j,k).isSingleValued()) {
            Print() << "norm = " << norm_arr(i,j,k,0) << "  " << norm_arr(i,j,k,1) << std::endl;
            if (   norm_arr(i,j,k,0) == eb_normal[0] 
                && norm_arr(i,j,k,1) == eb_normal[1]
#if (AMREX_SPACEDIM == 3)
                && norm_arr(i,j,k,2) == eb_normal[2]
#endif
               )
            {
               vel_arr(i,j,k,0) = eb_velocity[0];
               vel_arr(i,j,k,1) = eb_velocity[1];
#if (AMREX_SPACEDIM == 3)
               vel_arr(i,j,k,2) = eb_velocity[2];
#endif
            }
          }
        });
     }

    // We make sure to only fill "nghost" ghost cells so we don't accidentally 
    // over-write good ghost cell values with unfilled ghost cell values 
    IntVect ng_vect(AMREX_D_DECL(nghost,nghost,nghost));
    vel.EnforcePeriodicity(0,AMREX_SPACEDIM,ng_vect,gm.periodicity());
}
