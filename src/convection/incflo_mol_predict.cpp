#include <incflo.H>
#include <incflo_convection_K.H>
#include <incflo_MAC_bcs.H>

using namespace amrex;

void incflo::predict_vels_on_faces (int lev, MultiFab& u_mac, MultiFab& v_mac,
                                    MultiFab& w_mac, MultiFab const& vel)
{
#ifdef AMREX_USE_EB
    auto const& fact = this->EBFactory(lev);
    auto const& flags = fact.getMultiEBCellFlagFab();
    auto const& fcent = fact.getFaceCent();
    auto const& ccent = fact.getCentroid();
#endif

    Box const& domain = Geom(lev).Domain();
    Vector<BCRec> const& h_bcrec = get_velocity_bcrec();
    BCRec const* d_bcrec = get_velocity_bcrec_device_ptr();

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
                predict_vels_on_faces_eb(lev,bx,ubx,vbx,wbx,u,v,w,vcc,flagarr,fcx,fcy,fcz,ccc);
            }
            else
#endif
            {
                predict_vels_on_faces(lev,ubx,vbx,wbx,u,v,w,vcc);
            }

            incflo_set_mac_bcs(domain,ubx,vbx,wbx,u,v,w,vcc,h_bcrec,d_bcrec);
        }
    }
}

void incflo::predict_vels_on_faces (int lev, Box const& ubx, Box const& vbx, Box const& wbx,
                                    Array4<Real> const& u, Array4<Real> const& v,
                                    Array4<Real> const& w, Array4<Real const> const& vcc)
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

    if ((extdir_ilo and domain_ilo >= ubx.smallEnd(0)-1) or
        (extdir_ihi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(ubx, [vcc,extdir_ilo,extdir_ihi,domain_ilo,domain_ihi,u]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
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

    if ((extdir_jlo and domain_jlo >= vbx.smallEnd(1)-1) or
        (extdir_jhi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, [vcc,extdir_jlo,extdir_jhi,domain_jlo,domain_jhi,v]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
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

    if ((extdir_klo and domain_klo >= wbx.smallEnd(2)-1) or
        (extdir_khi and domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(wbx, [vcc,extdir_klo,extdir_khi,domain_klo,domain_khi,w]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
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
        });
    }
}
