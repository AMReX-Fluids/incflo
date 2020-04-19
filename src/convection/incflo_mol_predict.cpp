#include <MOL.H>
#include <incflo_slopes_K.H>
#include <incflo_MAC_bcs.H>

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir (BCRec const* bcrec, int ncomp, int dir)
    {
        std::pair<bool,bool> r{false,false};
        for (int n = 0; n < ncomp; ++n) {
            r.first = r.first or bcrec[n].lo(dir) == BCType::ext_dir;
            r.second = r.second or bcrec[n].hi(dir) == BCType::ext_dir;
        }
        return r;
    }
}

void 
mol::predict_vels_on_faces (int lev, MultiFab& u_mac, MultiFab& v_mac,
                            MultiFab& w_mac, MultiFab const& vel,
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

    Box const& domain = geom[lev].Domain();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(vel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& ubx = mfi.nodaltilebox(0);
            Box const& vbx = mfi.nodaltilebox(1);
            Box const& wbx = mfi.nodaltilebox(2);
            Array4<Real> const& u = u_mac.array(mfi);
            Array4<Real> const& v = v_mac.array(mfi);
            Array4<Real> const& w = w_mac.array(mfi);
            Array4<Real const> const& vcc = vel.const_array(mfi);
#ifdef AMREX_USE_EB
            Box const& bx = mfi.tilebox();
            EBCellFlagFab const& flagfab = flags[mfi];
            Array4<EBCellFlag const> const& flagarr = flagfab.const_array();
            auto const typ = flagfab.getType(amrex::grow(bx,1));
            if (typ == FabType::covered)
            {
                amrex::ParallelFor(ubx, vbx, wbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { u(i,j,k) = 0.0; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { v(i,j,k) = 0.0; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { w(i,j,k) = 0.0; });
            }
            else if (typ == FabType::singlevalued)
            {
                Array4<Real const> const& fcx = fcent[0]->const_array(mfi);
                Array4<Real const> const& fcy = fcent[1]->const_array(mfi);
                Array4<Real const> const& fcz = fcent[2]->const_array(mfi);
                Array4<Real const> const& ccc = ccent.const_array(mfi);
                predict_vels_on_faces_eb(lev,bx,ubx,vbx,wbx,
                                         u,v,w,vcc,flagarr,fcx,fcy,fcz,ccc,
                                         h_bcrec,d_bcrec,geom);
            }
            else
#endif
            {
                predict_vels_on_faces(lev,ubx,vbx,wbx,u,v,w,vcc,h_bcrec,d_bcrec,geom);
            }

            incflo_set_mac_bcs(domain,ubx,vbx,wbx,u,v,w,vcc,h_bcrec,d_bcrec);
        }
    }
}

void 
mol::predict_vels_on_faces (int lev, Box const& ubx, Box const& vbx, Box const& wbx,
                            Array4<Real> const& u, Array4<Real> const& v,
                            Array4<Real> const& w, Array4<Real const> const& vcc,
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
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir(h_bcrec.data(), ncomp, static_cast<int>(Direction::x));
    bool has_extdir_lo = extdir_lohi.first;
    bool has_extdir_hi = extdir_lohi.second;

    if ((has_extdir_lo and domain_ilo >= ubx.smallEnd(0)-1) or
        (has_extdir_hi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(ubx, [vcc,domain_ilo,domain_ihi,u,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_ilo = d_bcrec[0].lo(0) == BCType::ext_dir;
            bool extdir_ihi = d_bcrec[0].hi(0) == BCType::ext_dir;

            const Real vcc_pls = vcc(i,j,k,0);
            const Real vcc_mns = vcc(i-1,j,k,0);

            Real upls = vcc_pls - 0.5 * incflo_xslope_extdir
                (i,j,k,0,vcc, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);

            Real umns = vcc_mns + 0.5 * incflo_xslope_extdir
                (i-1,j,k,0,vcc, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);

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

            if (extdir_ilo and i == domain_ilo) {
                u_val = vcc_mns;
            } else if (extdir_ihi and i == domain_ihi+1) {
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

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    extdir_lohi = has_extdir(h_bcrec.data(), ncomp, static_cast<int>(Direction::y));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;

    if ((has_extdir_lo and domain_jlo >= vbx.smallEnd(1)-1) or
        (has_extdir_hi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, [vcc,domain_jlo,domain_jhi,v,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_jlo = d_bcrec[1].lo(1) == BCType::ext_dir;
            bool extdir_jhi = d_bcrec[1].hi(1) == BCType::ext_dir;

            const Real vcc_pls = vcc(i,j,k,1);
            const Real vcc_mns = vcc(i,j-1,k,1);

            Real vpls = vcc_pls - 0.5 * incflo_yslope_extdir
                (i,j,k,1,vcc, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
            Real vmns = vcc_mns + 0.5 * incflo_yslope_extdir
                (i,j-1,k,1,vcc, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);

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

            if (extdir_jlo and j == domain_jlo) {
                v_val = vcc_mns;
            } else if (extdir_jhi and j == domain_jhi+1) {
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

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    extdir_lohi = has_extdir(h_bcrec.data(), ncomp, static_cast<int>(Direction::z));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;

    if ((has_extdir_lo and domain_klo >= wbx.smallEnd(2)-1) or
        (has_extdir_hi and domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(wbx, [vcc,domain_klo,domain_khi,w,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_klo = d_bcrec[2].lo(2) == BCType::ext_dir;
            bool extdir_khi = d_bcrec[2].hi(2) == BCType::ext_dir;

            const Real vcc_pls = vcc(i,j,k,2);
            const Real vcc_mns = vcc(i,j,k-1,2);

            Real wpls = vcc_pls - 0.5 * incflo_zslope_extdir
                (i,j,k,2,vcc, extdir_klo, extdir_khi, domain_klo, domain_khi);
            Real wmns = vcc_mns + 0.5 * incflo_zslope_extdir(
                i,j,k-1,2,vcc, extdir_klo, extdir_khi, domain_klo, domain_khi);

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

            if (extdir_klo and k == domain_klo) {
                w_val = vcc_mns;
            } else if (extdir_khi and k == domain_khi+1) {
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
}
