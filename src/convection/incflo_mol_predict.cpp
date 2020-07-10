#include <MOL.H>
#include <incflo_slopes_K.H>

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

void 
mol::predict_vels_on_faces (int lev, 
                            AMREX_D_DECL(MultiFab& u_mac, 
                                         MultiFab& v_mac,
                                         MultiFab& w_mac), 
                            MultiFab const& vel,
                            Vector<BCRec> const& h_bcrec,
                                   BCRec  const* d_bcrec,
#ifdef AMREX_USE_EB
                            EBFArrayBoxFactory const* ebfact,
#endif
                            Vector<Geometry> geom)
{
#ifdef AMREX_USE_EB
    auto const& flags = ebfact->getMultiEBCellFlagFab();
    auto const& fcent = ebfact->getFaceCent();
    auto const& ccent = ebfact->getCentroid();
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(vel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            AMREX_D_TERM(Box const& ubx = mfi.nodaltilebox(0);,
                         Box const& vbx = mfi.nodaltilebox(1);,
                         Box const& wbx = mfi.nodaltilebox(2););
            AMREX_D_TERM(Array4<Real> const& u = u_mac.array(mfi);,
                         Array4<Real> const& v = v_mac.array(mfi);,
                         Array4<Real> const& w = w_mac.array(mfi););
            Array4<Real const> const& vcc = vel.const_array(mfi);
#ifdef AMREX_USE_EB
            Box const& bx = mfi.tilebox();
            EBCellFlagFab const& flagfab = flags[mfi];
            Array4<EBCellFlag const> const& flagarr = flagfab.const_array();
            auto const typ = flagfab.getType(amrex::grow(bx,2));
            if (typ == FabType::covered)
            {
#if (AMREX_SPACEDIM == 3)
                amrex::ParallelFor(ubx, vbx, wbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { u(i,j,k) = 0.0; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { v(i,j,k) = 0.0; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { w(i,j,k) = 0.0; });
#else
                amrex::ParallelFor(ubx, vbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { u(i,j,k) = 0.0; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { v(i,j,k) = 0.0; });
#endif
            }
            else if (typ == FabType::singlevalued)
            {
                AMREX_D_TERM(Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
                             Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
                             Array4<Real const> const& fcz = fcent[2]->const_array(mfi););
                Array4<Real const> const& ccc = ccent.const_array(mfi);
                predict_vels_on_faces_eb(lev,bx,AMREX_D_DECL(ubx,vbx,wbx),
                                         AMREX_D_DECL(u,v,w),vcc,flagarr,AMREX_D_DECL(fcx,fcy,fcz),ccc,
                                         h_bcrec,d_bcrec,geom);
            }
            else
#endif
            {
                predict_vels_on_faces(lev,AMREX_D_DECL(ubx,vbx,wbx),AMREX_D_DECL(u,v,w),vcc,h_bcrec,d_bcrec,geom);
            }
        }
    }
}

void 
mol::predict_vels_on_faces (int lev, 
                            AMREX_D_DECL(Box const& ubx, 
                                         Box const& vbx, 
                                         Box const& wbx),
                            AMREX_D_DECL(Array4<Real> const& u, 
                                         Array4<Real> const& v,
                                         Array4<Real> const& w), 
                            Array4<Real const> const& vcc,
                            Vector<BCRec> const& h_bcrec,
                                   BCRec  const* d_bcrec,
                            Vector<Geometry> geom)
{
    constexpr Real small_vel = 1.e-10;

    int ncomp = AMREX_SPACEDIM; // This is only used because h_bcrec and d_bcrec hold the 
                                // bc's for all three velocity components

    const Box& domain_box = geom[lev].Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);
#endif

    // At an ext_dir or hoextrap boundary, 
    //    the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo = extdir_lohi.first;
    bool has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_ilo >= ubx.smallEnd(0)-1) or
        (has_extdir_or_ho_hi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(ubx, [vcc,domain_ilo,domain_ihi,u,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_or_ho_ilo = (d_bcrec[0].lo(0) == BCType::ext_dir) or
                                    (d_bcrec[0].lo(0) == BCType::hoextrap);
            bool extdir_or_ho_ihi = (d_bcrec[0].hi(0) == BCType::ext_dir) or
                                    (d_bcrec[0].hi(0) == BCType::hoextrap);

            const Real vcc_pls = vcc(i,j,k,0);
            const Real vcc_mns = vcc(i-1,j,k,0);

            Real upls = vcc_pls - 0.5 * incflo_xslope_extdir
                (i,j,k,0,vcc, extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);

            Real umns = vcc_mns + 0.5 * incflo_xslope_extdir
                (i-1,j,k,0,vcc, extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);

            Real u_val(0);

            if (umns >= 0.0 or upls <= 0.0) {
                
                Real avg = 0.5 * (upls + umns);

                if (avg >= small_vel) {
                    u_val = umns;
                }
                else if (avg <= -small_vel){
                    u_val = upls;
                }
            }

            if (i == domain_ilo && (d_bcrec[0].lo(0) == BCType::ext_dir)) {
                u_val = vcc_mns;
            } else if (i == domain_ihi+1 && (d_bcrec[0].hi(0) == BCType::ext_dir)) {
                u_val = vcc_pls;
            }

            u(i,j,k) = u_val;
        });
    }
    else
    {
        amrex::ParallelFor(ubx, [vcc,u]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real upls = vcc(i  ,j,k,0) - 0.5 * incflo_xslope(i  ,j,k,0,vcc);
            Real umns = vcc(i-1,j,k,0) + 0.5 * incflo_xslope(i-1,j,k,0,vcc);

            Real u_val(0);

            if (umns >= 0.0 or upls <= 0.0) {
                
                Real avg = 0.5 * (upls + umns);
                
                if (avg >= small_vel) {
                    u_val = umns;
                }
                else if (avg <= -small_vel){
                    u_val = upls;
                }
            }

            u(i,j,k) = u_val;
        });
    }

    // At an ext_dir or hoextrap boundary, 
    //    the boundary value is on the face, not cell center.
    extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::y));
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_jlo >= vbx.smallEnd(1)-1) or
        (has_extdir_or_ho_hi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, [vcc,domain_jlo,domain_jhi,v,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_or_ho_jlo = (d_bcrec[1].lo(1) == BCType::ext_dir) or
                                    (d_bcrec[1].lo(1) == BCType::hoextrap);
            bool extdir_or_ho_jhi = (d_bcrec[1].hi(1) == BCType::ext_dir) or
                                    (d_bcrec[1].hi(1) == BCType::hoextrap);

            const Real vcc_pls = vcc(i,j,k,1);
            const Real vcc_mns = vcc(i,j-1,k,1);

            Real vpls = vcc_pls - 0.5 * incflo_yslope_extdir
                (i,j,k,1,vcc, extdir_or_ho_jlo, extdir_or_ho_jhi, domain_jlo, domain_jhi);
            Real vmns = vcc_mns + 0.5 * incflo_yslope_extdir
                (i,j-1,k,1,vcc, extdir_or_ho_jlo, extdir_or_ho_jhi, domain_jlo, domain_jhi);

            Real v_val(0);

            if (vmns >= 0.0 or vpls <= 0.0) {
                Real avg = 0.5 * (vpls + vmns);

                if (avg >= small_vel) {
                    v_val = vmns;
                }
                else if (avg <= -small_vel){
                    v_val = vpls;
                }
            }

            if (j == domain_jlo && (d_bcrec[1].lo(1) == BCType::ext_dir)) {
                v_val = vcc_mns;
            } else if (j == domain_jhi+1 && (d_bcrec[1].hi(1) == BCType::ext_dir)) {
                v_val = vcc_pls;
            }

            v(i,j,k) = v_val;
        });
    }
    else
    {
        amrex::ParallelFor(vbx, [vcc,v]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real vpls = vcc(i,j  ,k,1) - 0.5 * incflo_yslope(i,j  ,k,1,vcc);
            Real vmns = vcc(i,j-1,k,1) + 0.5 * incflo_yslope(i,j-1,k,1,vcc);

            Real v_val(0);

            if (vmns >= 0.0 or vpls <= 0.0) {
                Real avg = 0.5 * (vpls + vmns);

                if (avg >= small_vel) {
                    v_val = vmns;
                }
                else if (avg <= -small_vel) {
                    v_val = vpls;
                }
            }

            v(i,j,k) = v_val;
        });
    }

#if (AMREX_SPACEDIM == 3)
    // At an ext_dir or hoextrap boundary, 
    //    the boundary value is on the face, not cell center.
    extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::z));
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_klo >= wbx.smallEnd(2)-1) or
        (has_extdir_or_ho_hi and domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(wbx, [vcc,domain_klo,domain_khi,w,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_or_ho_klo = (d_bcrec[2].lo(2) == BCType::ext_dir) or
                                    (d_bcrec[2].lo(2) == BCType::hoextrap);
            bool extdir_or_ho_khi = (d_bcrec[2].hi(2) == BCType::ext_dir) or
                                    (d_bcrec[2].hi(2) == BCType::hoextrap);

            const Real vcc_pls = vcc(i,j,k,2);
            const Real vcc_mns = vcc(i,j,k-1,2);

            Real wpls = vcc_pls - 0.5 * incflo_zslope_extdir
                (i,j,k,2,vcc, extdir_or_ho_klo, extdir_or_ho_khi, domain_klo, domain_khi);
            Real wmns = vcc_mns + 0.5 * incflo_zslope_extdir(
                i,j,k-1,2,vcc, extdir_or_ho_klo, extdir_or_ho_khi, domain_klo, domain_khi);

            Real w_val(0);

            if (wmns >= 0.0 or wpls <= 0.0) {
                Real avg = 0.5 * (wpls + wmns);

                if (avg >= small_vel) {
                    w_val = wmns;
                }
                else if (avg <= -small_vel) {
                    w_val = wpls;
                }
            }

            if (k == domain_klo && (d_bcrec[2].lo(2) == BCType::ext_dir)) {
                w_val = vcc_mns;
            } else if (k == domain_khi+1 && (d_bcrec[2].hi(2) == BCType::ext_dir)) {
                w_val = vcc_pls;
            }

            w(i,j,k) = w_val;
        });
    }
    else
    {
        amrex::ParallelFor(wbx, [vcc,w]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real wpls = vcc(i,j,k  ,2) - 0.5 * incflo_zslope(i,j,k  ,2,vcc);
            Real wmns = vcc(i,j,k-1,2) + 0.5 * incflo_zslope(i,j,k-1,2,vcc);

            Real w_val(0);

            if (wmns >= 0.0 or wpls <= 0.0) {
                Real avg = 0.5 * (wpls + wmns);

                if (avg >= small_vel) {
                    w_val = wmns;
                }
                else if (avg <= -small_vel) {
                    w_val = wpls;
                }
            }

            w(i,j,k) = w_val;
        });
    }
#endif
}
