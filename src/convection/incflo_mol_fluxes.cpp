#include <incflo_slopes_K.H>
#include <MOL.H>
#include <utility>

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
mol::compute_convective_fluxes (int lev, Box const& bx, int ncomp,
                                AMREX_D_DECL(Array4<Real> const& fx,
                                             Array4<Real> const& fy,
                                             Array4<Real> const& fz),
                                Array4<Real const> const& q,
                                AMREX_D_DECL(Array4<Real const> const& umac,
                                             Array4<Real const> const& vmac,
                                             Array4<Real const> const& wmac),
                                BCRec const* h_bcrec, BCRec const* d_bcrec,
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

    AMREX_D_TERM(Box const& xbx = amrex::surroundingNodes(bx,0);,
                 Box const& ybx = amrex::surroundingNodes(bx,1);,
                 Box const& zbx = amrex::surroundingNodes(bx,2););

    // At an ext_dir or hoextrap boundary, 
    //    the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir_or_ho(h_bcrec, ncomp, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo = extdir_lohi.first;
    bool has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_ilo >= xbx.smallEnd(0)-1) or
        (has_extdir_or_ho_hi and domain_ihi <= xbx.bigEnd(0)))
    {
        amrex::ParallelFor(xbx, ncomp, [d_bcrec,q,domain_ilo,domain_ihi,umac,small_vel,fx]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            bool extdir_or_ho_ilo = (d_bcrec[n].lo(0) == BCType::ext_dir) or
                                    (d_bcrec[n].lo(0) == BCType::hoextrap);
            bool extdir_or_ho_ihi = (d_bcrec[n].hi(0) == BCType::ext_dir) or
                                    (d_bcrec[n].hi(0) == BCType::hoextrap);
            Real qs;
            if (i <= domain_ilo && (d_bcrec[n].lo(0) == BCType::ext_dir)) {
                qs = q(domain_ilo-1,j,k,n);
            } else if (i >= domain_ihi+1 && (d_bcrec[n].hi(0) == BCType::ext_dir)) {
                qs = q(domain_ihi+1,j,k,n);
            } else {
                Real qpls = q(i,j,k,n) - 0.5 * incflo_xslope_extdir
                    (i,j,k,n,q, extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);
                Real qmns = q(i-1,j,k,n) + 0.5 * incflo_xslope_extdir
                    (i-1,j,k,n,q, extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);
                if (umac(i,j,k) > small_vel) {
                    qs = qmns;
                } else if (umac(i,j,k) < -small_vel) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }
            }
            fx(i,j,k,n) = qs * umac(i,j,k);
        });
    }
    else
    {
        amrex::ParallelFor(xbx, ncomp, [q,umac,small_vel,fx]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls = q(i  ,j,k,n) - 0.5 * incflo_xslope(i  ,j,k,n,q);
            Real qmns = q(i-1,j,k,n) + 0.5 * incflo_xslope(i-1,j,k,n,q);
            Real qs;
            if (umac(i,j,k) > small_vel) {
                qs = qmns;
            } else if (umac(i,j,k) < -small_vel) {
                qs = qpls;
            } else {
                qs = 0.5*(qmns+qpls);
            }
            fx(i,j,k,n) = qs * umac(i,j,k);
        });
    }

    // At an ext_dir or hoextrap boundary, 
    //    the boundary value is on the face, not cell center.
    extdir_lohi = has_extdir_or_ho(h_bcrec, ncomp,  static_cast<int>(Direction::y));
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_jlo >= ybx.smallEnd(1)-1) or
        (has_extdir_or_ho_hi and domain_jhi <= ybx.bigEnd(1)))
    {
        amrex::ParallelFor(ybx, ncomp, [d_bcrec,q,domain_jlo,domain_jhi,vmac,small_vel,fy]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            bool extdir_or_ho_jlo = (d_bcrec[n].lo(1) == BCType::ext_dir) or
                                    (d_bcrec[n].lo(1) == BCType::hoextrap);
            bool extdir_or_ho_jhi = (d_bcrec[n].hi(1) == BCType::ext_dir) or
                                    (d_bcrec[n].hi(1) == BCType::hoextrap);
            Real qs;
            if (j <= domain_jlo && (d_bcrec[n].lo(1) == BCType::ext_dir)) {
                qs = q(i,domain_jlo-1,k,n);
            } else if (j >= domain_jhi+1 && (d_bcrec[n].hi(1) == BCType::ext_dir)) {
                qs = q(i,domain_jhi+1,k,n);
            } else {
                Real qpls = q(i,j,k,n) - 0.5 * incflo_yslope_extdir
                    (i,j,k,n,q, extdir_or_ho_jlo, extdir_or_ho_jhi, domain_jlo, domain_jhi);
                Real qmns = q(i,j-1,k,n) + 0.5 * incflo_yslope_extdir
                    (i,j-1,k,n,q, extdir_or_ho_jlo, extdir_or_ho_jhi, domain_jlo, domain_jhi);
                if (vmac(i,j,k) > small_vel) {
                    qs = qmns;
                } else if (vmac(i,j,k) < -small_vel) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }
            }
            fy(i,j,k,n) = qs * vmac(i,j,k);
        });
    }
    else
    {
        amrex::ParallelFor(ybx, ncomp, [q,vmac,small_vel,fy]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls = q(i,j  ,k,n) - 0.5 * incflo_yslope(i,j  ,k,n,q);
            Real qmns = q(i,j-1,k,n) + 0.5 * incflo_yslope(i,j-1,k,n,q);
            Real qs;
            if (vmac(i,j,k) > small_vel) {
                qs = qmns;
            } else if (vmac(i,j,k) < -small_vel) {
                qs = qpls;
            } else {
                qs = 0.5*(qmns+qpls);
            }
            fy(i,j,k,n) = qs * vmac(i,j,k);
        });
    }

#if (AMREX_SPACEDIM == 3)
    // At an ext_dir or hoextrap boundary, 
    //    the boundary value is on the face, not cell center.
    extdir_lohi = has_extdir_or_ho(h_bcrec, ncomp, static_cast<int>(Direction::z));
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_klo >= zbx.smallEnd(2)-1) or
        (has_extdir_or_ho_hi and domain_khi <= zbx.bigEnd(2)))
    {
        amrex::ParallelFor(zbx, ncomp, [d_bcrec,q,domain_klo,domain_khi,wmac,small_vel,fz]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            bool extdir_or_ho_klo =   (d_bcrec[n].lo(2) == BCType::ext_dir) or
                                      (d_bcrec[n].lo(2) == BCType::hoextrap);
            bool extdir_or_ho_khi =   (d_bcrec[n].hi(2) == BCType::ext_dir) or
                                      (d_bcrec[n].hi(2) == BCType::hoextrap);
            Real qs;
            if (k <= domain_klo && (d_bcrec[n].lo(2) == BCType::ext_dir)) {
                qs = q(i,j,domain_klo-1,n);
            } else if (k >= domain_khi+1 && (d_bcrec[n].hi(2) == BCType::ext_dir)) {
                qs = q(i,j,domain_khi+1,n);
            } else {
                Real qpls = q(i,j,k,n) - 0.5 * incflo_zslope_extdir
                    (i,j,k,n,q, extdir_or_ho_klo, extdir_or_ho_khi, domain_klo, domain_khi);
                Real qmns = q(i,j,k-1,n) + 0.5 * incflo_zslope_extdir(
                    i,j,k-1,n,q, extdir_or_ho_klo, extdir_or_ho_khi, domain_klo, domain_khi);
                if (wmac(i,j,k) > small_vel) {
                    qs = qmns;
                } else if (wmac(i,j,k) < -small_vel) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }
            }
            fz(i,j,k,n) = qs * wmac(i,j,k);
        });
    }
    else
    {
        amrex::ParallelFor(zbx, ncomp, [q,wmac,small_vel,fz]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls = q(i,j,k  ,n) - 0.5 * incflo_zslope(i,j,k  ,n,q);
            Real qmns = q(i,j,k-1,n) + 0.5 * incflo_zslope(i,j,k-1,n,q);
            Real qs;
            if (wmac(i,j,k) > small_vel) {
                qs = qmns;
            } else if (wmac(i,j,k) < -small_vel) {
                qs = qpls;
            } else {
                qs = 0.5*(qmns+qpls);
            }
            fz(i,j,k,n) = qs * wmac(i,j,k);
        });
    }
#endif
}

