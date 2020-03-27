#include <incflo_convection_K.H>
#include <incflo.H>

using namespace amrex;

#ifdef AMREX_USE_EB
void incflo::predict_vels_on_faces_eb (int lev, Box const& ccbx,
                                       Box const& ubx, Box const& vbx, Box const& wbx,
                                       Array4<Real> const& u, Array4<Real> const& v,
                                       Array4<Real> const& w, Array4<Real const> const& vcc,
                                       Array4<EBCellFlag const> const& flag,
                                       Array4<Real const> const& fcx,
                                       Array4<Real const> const& fcy,
                                       Array4<Real const> const& fcz,
                                       Array4<Real const> const& ccc)
{
    constexpr Real small_vel = 1.e-10;

    const Box& domain_box = geom[lev].Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);

    auto const bc_ilo = m_bc_type[Orientation(Direction::x,Orientation::low)];
    auto const bc_ihi = m_bc_type[Orientation(Direction::x,Orientation::high)];
    auto const bc_jlo = m_bc_type[Orientation(Direction::y,Orientation::low)];
    auto const bc_jhi = m_bc_type[Orientation(Direction::y,Orientation::high)];
    auto const bc_klo = m_bc_type[Orientation(Direction::z,Orientation::low)];
    auto const bc_khi = m_bc_type[Orientation(Direction::z,Orientation::high)];

    bool extdir_ilo = (bc_ilo == BC::mass_inflow) or (bc_ilo == BC::no_slip_wall);
    bool extdir_ihi = (bc_ihi == BC::mass_inflow) or (bc_ihi == BC::no_slip_wall);
    bool extdir_jlo = (bc_jlo == BC::mass_inflow) or (bc_jlo == BC::no_slip_wall);
    bool extdir_jhi = (bc_jhi == BC::mass_inflow) or (bc_jhi == BC::no_slip_wall);
    bool extdir_klo = (bc_klo == BC::mass_inflow) or (bc_klo == BC::no_slip_wall);
    bool extdir_khi = (bc_khi == BC::mass_inflow) or (bc_khi == BC::no_slip_wall);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.

    // ****************************************************************************
    // Predict to x-faces
    // ****************************************************************************
    if ((extdir_ilo and domain_ilo >= ubx.smallEnd(0)-1) or
        (extdir_ihi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(Box(ubx),
        [u,vcc,flag,fcx,ccc,extdir_ilo,extdir_ihi,extdir_jlo,extdir_jhi,extdir_klo,extdir_khi,
         domain_ilo,domain_ihi,domain_jlo,domain_jhi,domain_klo,domain_khi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real u_val(0);

            if (flag(i,j,k).isConnected(-1,0,0))
            {
               Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
               Real zf = fcx(i,j,k,1);

               Real delta_x = 0.5 + ccc(i,j,k,0);
               Real delta_y = yf  - ccc(i,j,k,1);
               Real delta_z = zf  - ccc(i,j,k,2);

               Real vcc_mns = vcc(i-1,j,k,0);
               Real vcc_pls = vcc(i,j,k,0);

               Real cc_umax = amrex::max(vcc_pls, vcc_mns);
               Real cc_umin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,0,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];

               upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

               delta_x = 0.5 - ccc(i-1,j,k,0);
               delta_y = yf  - ccc(i-1,j,k,1);
               delta_z = zf  - ccc(i-1,j,k,2);

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i-1,j,k,0,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];

               umns = amrex::max(amrex::min(umns, cc_umax), cc_umin);

               if ( umns >= 0.0 or upls <= 0.0 ) {
                  Real avg = 0.5 * ( upls + umns );

                  if (avg >= small_vel) {
                    u_val = umns;
                  }
                  else if (avg <= -small_vel) {
                    u_val = upls;
                  }
               }

               if (extdir_ilo and i == domain_ilo) {
                   u_val = vcc_mns;
               } else if (extdir_ihi and i == domain_ihi+1) {
                   u_val = vcc_pls;
               }
            }
            
            u(i,j,k) = u_val;
        });
    }
    else
    {
        amrex::ParallelFor(Box(ubx),
        [u,vcc,flag,fcx,ccc]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real u_val(0);

            if (flag(i,j,k).isConnected(-1,0,0))
            {
               Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
               Real zf = fcx(i,j,k,1);

               Real delta_x = 0.5 + ccc(i,j,k,0);
               Real delta_y = yf  - ccc(i,j,k,1);
               Real delta_z = zf  - ccc(i,j,k,2);

               const Real vcc_mns = vcc(i-1,j,k,0);
               const Real vcc_pls = vcc(i,j,k,0);

               Real cc_umax = amrex::max(vcc_pls, vcc_mns);
               Real cc_umin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "0" of vcc
               const auto slopes_eb_hi = incflo_slopes_eb(i,j,k,0,vcc,ccc,flag);

               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];

               upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

               delta_x = 0.5 - ccc(i-1,j,k,0);
               delta_y = yf  - ccc(i-1,j,k,1);
               delta_z = zf  - ccc(i-1,j,k,2);

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_lo = incflo_slopes_eb(i-1,j,k,0,vcc,ccc,flag);

               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];

               umns = amrex::max(amrex::min(umns, cc_umax), cc_umin);

               if ( umns >= 0.0 or upls <= 0.0 ) {
                  Real avg = 0.5 * ( upls + umns );

                  if (avg >= small_vel) {
                    u_val = umns;
                  }
                  else if (avg <= -small_vel) {
                    u_val = upls;
                  }
               }
            }

            u(i,j,k) = u_val;
        });
    }

    // ****************************************************************************
    // Predict to y-faces
    // ****************************************************************************
    if ((extdir_jlo and domain_jlo >= vbx.smallEnd(1)-1) or
        (extdir_jhi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(Box(vbx),
        [v,vcc,flag,fcy,ccc,extdir_ilo,extdir_ihi,extdir_jlo,extdir_jhi,extdir_klo,extdir_khi,
         domain_ilo,domain_ihi,domain_jlo,domain_jhi,domain_klo,domain_khi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real v_val(0);

            if (flag(i,j,k).isConnected(0,-1,0))
            {
               Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
               Real zf = fcy(i,j,k,1);

               Real delta_x = xf  - ccc(i,j,k,0);
               Real delta_y = 0.5 + ccc(i,j,k,1);
               Real delta_z = zf  - ccc(i,j,k,2);

               const Real vcc_mns = vcc(i,j-1,k,1);
               const Real vcc_pls = vcc(i,j,k,1);

               Real cc_vmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_vmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,1,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];

               vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

               delta_x = xf  - ccc(i,j-1,k,0);
               delta_y = 0.5 - ccc(i,j-1,k,1);
               delta_z = zf  - ccc(i,j-1,k,2);

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i,j-1,k,1,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];

               vmns = amrex::max(amrex::min(vmns, cc_vmax), cc_vmin);

               if ( vmns >= 0.0 or vpls <= 0.0 ) {
                  Real avg = 0.5 * ( vpls + vmns );

                  if (avg >= small_vel) {
                    v_val = vmns;
                  }
                  else if (avg <= -small_vel) {
                    v_val = vpls;
                  }
               }

               if (extdir_jlo and j == domain_jlo) {
                   v_val = vcc_mns;
               } 
               else if (extdir_jhi and j == domain_jhi+1) {
                   v_val = vcc_pls;
               }
            }

            v(i,j,k) = v_val;
        });
    }
    else
    {
        amrex::ParallelFor(Box(vbx),
        [v,vcc,flag,fcy,ccc] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real v_val(0);

            if (flag(i,j,k).isConnected(0,-1,0))
            {
               Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
               Real zf = fcy(i,j,k,1);

               Real delta_x = xf  - ccc(i,j,k,0);
               Real delta_y = 0.5 + ccc(i,j,k,1);
               Real delta_z = zf  - ccc(i,j,k,2);

               const Real vcc_mns = vcc(i,j-1,k,1);
               const Real vcc_pls = vcc(i,j,k,1);

               Real cc_vmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_vmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "1" of vcc
               const auto slopes_eb_hi = incflo_slopes_eb(i,j,k,1,vcc,ccc,flag);

               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];

               vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

               delta_x = xf  - ccc(i,j-1,k,0);
               delta_y = 0.5 - ccc(i,j-1,k,1);
               delta_z = zf  - ccc(i,j-1,k,2);

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_lo = incflo_slopes_eb(i,j-1,k,1,vcc,ccc,flag);

               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
                                          
               vmns = amrex::max(amrex::min(vmns, cc_vmax), cc_vmin);

               if ( vmns >= 0.0 or vpls <= 0.0 ) {
                  Real avg = 0.5 * ( vpls + vmns );

                  if (avg >= small_vel) {
                    v(i,j,k) = vmns;
                  }
                  else if (avg <= -small_vel) {
                    v(i,j,k) = vpls;
                  }
               }
            }

            v(i,j,k) = v_val;
        });
    }

    // ****************************************************************************
    // Predict to z-faces
    // ****************************************************************************
    if ((extdir_klo and domain_klo >= wbx.smallEnd(2)-1) or
        (extdir_khi and domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(Box(wbx),
        [w,vcc,flag,fcz,ccc,extdir_ilo,extdir_ihi,extdir_jlo,extdir_jhi,extdir_klo,extdir_khi,
         domain_ilo,domain_ihi,domain_jlo,domain_jhi,domain_klo,domain_khi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real w_val(0);

            if (flag(i,j,k).isConnected(0,0,-1))
            {
               Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
               Real yf = fcz(i,j,k,1);

               Real delta_x = xf  - ccc(i,j,k,0);
               Real delta_y = yf  - ccc(i,j,k,1);
               Real delta_z = 0.5 + ccc(i,j,k,2);

               const Real vcc_mns = vcc(i,j,k-1,2);
               const Real vcc_pls = vcc(i,j,k,2);

               Real cc_wmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_wmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,2,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real wpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   - delta_z * slopes_eb_hi[2];

               wpls = amrex::max(amrex::min(wpls, cc_wmax), cc_wmin);

               delta_x = xf  - ccc(i,j,k-1,0);
               delta_y = yf  - ccc(i,j,k-1,1);
               delta_z = 0.5 - ccc(i,j,k-1,2);

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i,j,k-1,2,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real wmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];

               wmns = amrex::max(amrex::min(wmns, cc_wmax), cc_wmin);

               if ( wmns >= 0.0 or wpls <= 0.0 ) {
                  Real avg = 0.5 * ( wpls + wmns );

                  if (avg >= small_vel) {
                    w_val = wmns;
                  }
                  else if (avg <= -small_vel) {
                    w_val = wpls;
                  }
               }

                if (extdir_klo and k == domain_klo) {
                    w_val = vcc_mns;
                }
                else if (extdir_khi and k == domain_khi+1) {
                    w_val = vcc_pls;
                }
            }

            w(i,j,k) = w_val;
        });
    }
    else
    {
        amrex::ParallelFor(Box(wbx),
        [w,vcc,flag,fcz,ccc] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real w_val(0);

            if (flag(i,j,k).isConnected(0,0,-1))
            {
               Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
               Real yf = fcz(i,j,k,1);

               Real delta_x = xf  - ccc(i,j,k,0);
               Real delta_y = yf  - ccc(i,j,k,1);
               Real delta_z = 0.5 + ccc(i,j,k,2);

               const Real vcc_mns = vcc(i,j,k-1,2);
               const Real vcc_pls = vcc(i,j,k,2);

               Real cc_wmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_wmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "2" of vcc
               const auto slopes_eb_hi = incflo_slopes_eb(i,j,k,2,vcc,ccc,flag);

               Real wpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   - delta_z * slopes_eb_hi[2];

               wpls = amrex::max(amrex::min(wpls, cc_wmax), cc_wmin);

               delta_x = xf  - ccc(i,j,k-1,0);
               delta_y = yf  - ccc(i,j,k-1,1);
               delta_z = 0.5 - ccc(i,j,k-1,2);

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_lo = incflo_slopes_eb(i,j,k-1,2,vcc,ccc,flag);

               Real wmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];

               wmns = amrex::max(amrex::min(wmns, cc_wmax), cc_wmin);

               if ( wmns >= 0.0 or wpls <= 0.0 ) {
                  Real avg = 0.5 * ( wpls + wmns );

                  if (avg >= small_vel) {
                    w_val = wmns;
                  }
                  else if (avg <= -small_vel) {
                    w_val = wpls;
                  }
               }
            }

            w(i,j,k) = w_val;
        });
    }
}
#endif
