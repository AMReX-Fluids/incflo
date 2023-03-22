#include <incflo.H>

using namespace amrex;

#ifdef AMREX_USE_EB
void
incflo::incflo_correct_small_cells (Vector<MultiFab*      > const& vel_in,
                                    AMREX_D_DECL(Vector<MultiFab const*> const& u_mac,
                                                 Vector<MultiFab const*> const& v_mac,
                                                 Vector<MultiFab const*> const& w_mac))
{
    BL_PROFILE("incflo::incflo_correct_small_cells");

    for (int lev = 0; lev <= finest_level; lev++)
    {
       for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
          // Tilebox
          const Box bx = mfi.tilebox();

          EBCellFlagFab const& flags = EBFactory(lev).getMultiEBCellFlagFab()[mfi];

          // Face-centered velocity components
          AMREX_D_TERM(const auto& umac_fab = (u_mac[lev])->array(mfi);,
                       const auto& vmac_fab = (v_mac[lev])->array(mfi);,
                       const auto& wmac_fab = (w_mac[lev])->array(mfi););

          if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
          {
            // do nothing
          }

          // No cut cells in this FAB
          else if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
          {
            // do nothing
          }

          // Cut cells in this FAB
          else
          {
             // Face-centered areas
             AMREX_D_TERM(const auto& apx_fab   = EBFactory(lev).getAreaFrac()[0]->const_array(mfi);,
                          const auto& apy_fab   = EBFactory(lev).getAreaFrac()[1]->const_array(mfi);,
                          const auto& apz_fab   = EBFactory(lev).getAreaFrac()[2]->const_array(mfi););

             const auto& vfrac_fab = EBFactory(lev).getVolFrac().const_array(mfi);

             const auto& ccvel_fab = vel_in[lev]->array(mfi);

             if (!m_eb_flow.enabled) {
                // This FAB has cut cells -- we define the centroid value in terms of the MAC velocities onfaces
                amrex::ParallelFor(bx,
                  [vfrac_fab,AMREX_D_DECL(apx_fab,apy_fab,apz_fab),ccvel_fab,AMREX_D_DECL(umac_fab,vmac_fab,wmac_fab)]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (vfrac_fab(i,j,k) > 0.0 && vfrac_fab(i,j,k) < 1.e-4)
                    {
                       AMREX_D_TERM(Real u_avg = (apx_fab(i,j,k) * umac_fab(i,j,k) + apx_fab(i+1,j,k) * umac_fab(i+1,j,k))
                                               / (apx_fab(i,j,k) + apx_fab(i+1,j,k));,
                                    Real v_avg = (apy_fab(i,j,k) * vmac_fab(i,j,k) + apy_fab(i,j+1,k) * vmac_fab(i,j+1,k))
                                               / (apy_fab(i,j,k) + apy_fab(i,j+1,k));,
                                    Real w_avg = (apz_fab(i,j,k) * wmac_fab(i,j,k) + apz_fab(i,j,k+1) * wmac_fab(i,j,k+1))
                                               / (apz_fab(i,j,k) + apz_fab(i,j,k+1)););

                       AMREX_D_TERM(ccvel_fab(i,j,k,0) = u_avg;,
                                    ccvel_fab(i,j,k,1) = v_avg;,
                                    ccvel_fab(i,j,k,2) = w_avg;);

                    }
                });
             } else { // EB has flow
                Array4<Real const> const& eb_vel = get_velocity_eb()[lev]->const_array(mfi);

                // This FAB has cut cells -- we define the centroid value in terms of the MAC velocities onfaces
                amrex::ParallelFor(bx,
                  [vfrac_fab,AMREX_D_DECL(apx_fab,apy_fab,apz_fab),ccvel_fab,AMREX_D_DECL(umac_fab,vmac_fab,wmac_fab),eb_vel]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (vfrac_fab(i,j,k) > 0.0 && vfrac_fab(i,j,k) < 1.e-4)
                    {
                       AMREX_D_TERM(const Real ucc = ccvel_fab(i,j,k,0);,
                                    const Real vcc = ccvel_fab(i,j,k,1);,
                                    const Real wcc = ccvel_fab(i,j,k,2););

                       AMREX_D_TERM(Real u_avg = (apx_fab(i,j,k) * umac_fab(i,j,k) + apx_fab(i+1,j,k) * umac_fab(i+1,j,k))
                                               / (apx_fab(i,j,k) + apx_fab(i+1,j,k));,
                                    Real v_avg = (apy_fab(i,j,k) * vmac_fab(i,j,k) + apy_fab(i,j+1,k) * vmac_fab(i,j+1,k))
                                               / (apy_fab(i,j,k) + apy_fab(i,j+1,k));,
                                    Real w_avg = (apz_fab(i,j,k) * wmac_fab(i,j,k) + apz_fab(i,j,k+1) * wmac_fab(i,j,k+1))
                                               / (apz_fab(i,j,k) + apz_fab(i,j,k+1)););


                       AMREX_D_TERM(ccvel_fab(i,j,k,0) = u_avg;,
                                    ccvel_fab(i,j,k,1) = v_avg;,
                                    ccvel_fab(i,j,k,2) = w_avg;);

                       AMREX_D_TERM(const Real u_eb = eb_vel(i,j,k,0);,
                                    const Real v_eb = eb_vel(i,j,k,1);,
                                    const Real w_eb = eb_vel(i,j,k,2););
#if (AMREX_SPACEDIM == 2)
                       const Real eb_vel_mag = std::sqrt(u_eb*u_eb + v_eb*v_eb);
#elif (AMREX_SPACEDIM == 3)
                       const Real eb_vel_mag = std::sqrt(u_eb*u_eb + v_eb*v_eb + w_eb*w_eb);
#endif
                       constexpr amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

                       if (eb_vel_mag > tolerance) {
                          if ( u_eb > tolerance ) {
                              ccvel_fab(i,j,k,0) = amrex::min(ucc, u_eb);
                          } else if( u_eb < -tolerance) {
                              ccvel_fab(i,j,k,0) = amrex::max(ucc, u_eb);
                          }
                          if ( v_eb > tolerance ) {
                              ccvel_fab(i,j,k,1) = amrex::min(vcc, v_eb);
                          } else if( v_eb < -tolerance) {
                              ccvel_fab(i,j,k,1) = amrex::max(vcc, v_eb);
                          }
#if (AMREX_SPACEDIM == 3)
                          if ( w_eb > tolerance ) {
                              ccvel_fab(i,j,k,2) = amrex::min(wcc, w_eb);
                          } else if( w_eb < -tolerance) {
                              ccvel_fab(i,j,k,2) = amrex::max(wcc, w_eb);
                          }
#endif
                       }
                    }
                });
             }
          }
       }
    }
}
#endif
