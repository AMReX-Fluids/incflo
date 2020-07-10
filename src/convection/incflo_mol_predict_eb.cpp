#include <incflo_slopes_K.H>
#include <MOL.H>

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir_or_ho (BCRec const* bcrec, int ncomp, int dir)
    {
        std::pair<bool,bool> r{false,false};
        for (int n = 0; n < ncomp; ++n) {
            r.first = r.first 
                 or (bcrec[n].lo(dir) == BCType::ext_dir)
                 or (bcrec[n].lo(dir) == BCType::hoextrap);
            r.second = r.second 
                 or (bcrec[n].hi(dir) == BCType::ext_dir)
                 or (bcrec[n].hi(dir) == BCType::hoextrap);
        }
        return r;
    }
}

#ifdef AMREX_USE_EB
void 
mol::predict_vels_on_faces_eb (int lev, Box const& ccbx,
                               AMREX_D_DECL(Box const& ubx, 
                                            Box const& vbx, 
                                            Box const& wbx),
                               AMREX_D_DECL(Array4<Real> const& u, 
                                            Array4<Real> const& v,
                                            Array4<Real> const& w), 
                               Array4<Real const> const& vcc,
                               Array4<EBCellFlag const> const& flag,
                               AMREX_D_DECL(Array4<Real const> const& fcx,
                                            Array4<Real const> const& fcy,
                                            Array4<Real const> const& fcz),
                               Array4<Real const> const& ccc,
                               Vector<BCRec> const& h_bcrec,
                                      BCRec  const* d_bcrec,
                               Vector<Geometry> geom)
{
    constexpr Real small_vel = 1.e-10;

    const Box& domain_box = geom[lev].Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);
#endif

    int ncomp = AMREX_SPACEDIM; // This is only used because h_bcrec and d_bcrec hold the
                                // bc's for all three velocity components

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo = extdir_lohi.first;
    bool has_extdir_or_ho_hi = extdir_lohi.second;

    // ****************************************************************************
    // Predict to x-faces
    // ****************************************************************************
    if ((has_extdir_or_ho_lo and domain_ilo >= ubx.smallEnd(0)-1) or
        (has_extdir_or_ho_hi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(Box(ubx),
        [u,vcc,flag,fcx,ccc,d_bcrec,
         AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
         AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi)]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real u_val(0);

            AMREX_D_TERM(bool extdir_or_ho_ilo = (d_bcrec[0].lo(0) == BCType::ext_dir) or
                                                 (d_bcrec[0].lo(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jlo = (d_bcrec[0].lo(1) == BCType::ext_dir) or
                                                 (d_bcrec[0].lo(1) == BCType::hoextrap);,
                         bool extdir_or_ho_klo = (d_bcrec[0].lo(2) == BCType::ext_dir) or
                                                 (d_bcrec[0].lo(2) == BCType::hoextrap););

            AMREX_D_TERM(bool extdir_or_ho_ihi = (d_bcrec[0].hi(0) == BCType::ext_dir) or
                                                 (d_bcrec[0].hi(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jhi = (d_bcrec[0].hi(1) == BCType::ext_dir) or
                                                 (d_bcrec[0].hi(1) == BCType::hoextrap);,
                         bool extdir_or_ho_khi = (d_bcrec[0].hi(2) == BCType::ext_dir) or
                                                 (d_bcrec[0].hi(2) == BCType::hoextrap););

            if (flag(i,j,k).isConnected(-1,0,0))
            {
               Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcx(i,j,k,1);
#endif
               AMREX_D_TERM(Real delta_x = 0.5 + ccc(i,j,k,0);,
                            Real delta_y = yf  - ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               Real vcc_mns = vcc(i-1,j,k,0);
               Real vcc_pls = vcc(i,j,k,0);

               Real cc_umax = amrex::max(vcc_pls, vcc_mns);
               Real cc_umin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,0,vcc,ccc,flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

#if (AMREX_SPACEDIM == 3)
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1];
#endif
               upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

               AMREX_D_TERM(delta_x = 0.5 - ccc(i-1,j,k,0);,
                            delta_y = yf  - ccc(i-1,j,k,1);,
                            delta_z = zf  - ccc(i-1,j,k,2););

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i-1,j,k,0,vcc,ccc,flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

#if (AMREX_SPACEDIM == 3)
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif
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

               if (i == domain_ilo && (d_bcrec[0].lo(0) == BCType::ext_dir)) {
                   u_val = vcc_mns;
               } else if (i == domain_ihi+1 && (d_bcrec[0].hi(0) == BCType::ext_dir)) {
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
#if (AMREX_SPACEDIM == 3)
               Real zf = fcx(i,j,k,1);
#endif

               AMREX_D_TERM(Real delta_x = 0.5 + ccc(i,j,k,0);,
                            Real delta_y = yf  - ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               const Real vcc_mns = vcc(i-1,j,k,0);
               const Real vcc_pls = vcc(i,j,k,0);

               Real cc_umax = amrex::max(vcc_pls, vcc_mns);
               Real cc_umin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "0" of vcc
               const auto slopes_eb_hi = incflo_slopes_eb(i,j,k,0,vcc,ccc,flag);

#if (AMREX_SPACEDIM == 3)
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1];
#endif
               upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

               AMREX_D_TERM(delta_x = 0.5 - ccc(i-1,j,k,0);,
                            delta_y = yf  - ccc(i-1,j,k,1);,
                            delta_z = zf  - ccc(i-1,j,k,2););

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_lo = incflo_slopes_eb(i-1,j,k,0,vcc,ccc,flag);

#if (AMREX_SPACEDIM == 3)
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif
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
    extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::y));
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_jlo >= vbx.smallEnd(1)-1) or
        (has_extdir_or_ho_hi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(Box(vbx),
        [v,vcc,flag,fcy,ccc,d_bcrec,
         AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
         AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi)]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real v_val(0);

            AMREX_D_TERM(bool extdir_or_ho_ilo = (d_bcrec[1].lo(0) == BCType::ext_dir) or
                                                 (d_bcrec[1].lo(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jlo = (d_bcrec[1].lo(1) == BCType::ext_dir) or
                                                 (d_bcrec[1].lo(1) == BCType::hoextrap);,
                         bool extdir_or_ho_klo = (d_bcrec[1].lo(2) == BCType::ext_dir) or
                                                 (d_bcrec[1].lo(2) == BCType::hoextrap););

            AMREX_D_TERM(bool extdir_or_ho_ihi = (d_bcrec[1].hi(0) == BCType::ext_dir) or
                                                 (d_bcrec[1].hi(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jhi = (d_bcrec[1].hi(1) == BCType::ext_dir) or
                                                 (d_bcrec[1].hi(1) == BCType::hoextrap);,
                         bool extdir_or_ho_khi = (d_bcrec[1].hi(2) == BCType::ext_dir) or
                                                 (d_bcrec[1].hi(2) == BCType::hoextrap););

            if (flag(i,j,k).isConnected(0,-1,0))
            {
               Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcy(i,j,k,1);
#endif

               AMREX_D_TERM(Real delta_x = xf  - ccc(i,j,k,0);,
                            Real delta_y = 0.5 + ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               const Real vcc_mns = vcc(i,j-1,k,1);
               const Real vcc_pls = vcc(i,j,k,1);

               Real cc_vmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_vmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,1,vcc,ccc,flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

#if (AMREX_SPACEDIM == 3)
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1];
#endif

               vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

               AMREX_D_TERM(delta_x = xf  - ccc(i,j-1,k,0);,
                            delta_y = 0.5 - ccc(i,j-1,k,1);,
                            delta_z = zf  - ccc(i,j-1,k,2););

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i,j-1,k,1,vcc,ccc,flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

#if (AMREX_SPACEDIM == 3)
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif

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

               if (j == domain_jlo && (d_bcrec[1].lo(1) == BCType::ext_dir)) {
                   v_val = vcc_mns;
               } 
               else if (j == domain_jhi+1 && (d_bcrec[1].hi(1) == BCType::ext_dir)) {
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
#if (AMREX_SPACEDIM == 3)
               Real zf = fcy(i,j,k,1);
#endif

               AMREX_D_TERM(Real delta_x = xf  - ccc(i,j,k,0);,
                            Real delta_y = 0.5 + ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               const Real vcc_mns = vcc(i,j-1,k,1);
               const Real vcc_pls = vcc(i,j,k,1);

               Real cc_vmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_vmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "1" of vcc
               const auto slopes_eb_hi = incflo_slopes_eb(i,j,k,1,vcc,ccc,flag);

#if (AMREX_SPACEDIM == 3)
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1];
#endif

               vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

               AMREX_D_TERM(delta_x = xf  - ccc(i,j-1,k,0);,
                            delta_y = 0.5 - ccc(i,j-1,k,1);,
                            delta_z = zf  - ccc(i,j-1,k,2););

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_lo = incflo_slopes_eb(i,j-1,k,1,vcc,ccc,flag);

#if (AMREX_SPACEDIM == 3)
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif
                                          
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
            }

            v(i,j,k) = v_val;
        });
    }

#if (AMREX_SPACEDIM == 3)
    // ****************************************************************************
    // Predict to z-faces
    // ****************************************************************************
    extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::z));
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_klo >= wbx.smallEnd(2)-1) or
        (has_extdir_or_ho_hi and domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(Box(wbx),
        [w,vcc,flag,fcz,ccc,d_bcrec,
         domain_ilo,domain_ihi,domain_jlo,domain_jhi,domain_klo,domain_khi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real w_val(0);

            bool extdir_or_ho_ilo = (d_bcrec[2].lo(0) == BCType::ext_dir) or
                                    (d_bcrec[2].lo(0) == BCType::hoextrap);
            bool extdir_or_ho_ihi = (d_bcrec[2].hi(0) == BCType::ext_dir) or
                                    (d_bcrec[2].hi(0) == BCType::hoextrap);

            bool extdir_or_ho_jlo = (d_bcrec[2].lo(1) == BCType::ext_dir) or
                                    (d_bcrec[2].lo(1) == BCType::hoextrap);
            bool extdir_or_ho_jhi = (d_bcrec[2].hi(1) == BCType::ext_dir) or
                                    (d_bcrec[2].hi(1) == BCType::hoextrap);

            bool extdir_or_ho_klo = (d_bcrec[2].lo(2) == BCType::ext_dir) or
                                    (d_bcrec[2].lo(2) == BCType::hoextrap);
            bool extdir_or_ho_khi = (d_bcrec[2].hi(2) == BCType::ext_dir) or
                                    (d_bcrec[2].hi(2) == BCType::hoextrap);

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
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

               Real wpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   - delta_z * slopes_eb_hi[2];

               wpls = amrex::max(amrex::min(wpls, cc_wmax), cc_wmin);

               delta_x = xf  - ccc(i,j,k-1,0);
               delta_y = yf  - ccc(i,j,k-1,1);
               delta_z = 0.5 - ccc(i,j,k-1,2);

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i,j,k-1,2,vcc,ccc,flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

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

                if (k == domain_klo && (d_bcrec[2].lo(2) == BCType::ext_dir)) {
                    w_val = vcc_mns;
                }
                else if (k == domain_khi+1 && (d_bcrec[2].hi(2) == BCType::ext_dir)) {
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
#endif
}
#endif
