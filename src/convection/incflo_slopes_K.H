#ifndef INCFLO_CONVECTION_K_H_
#define INCFLO_CONVECTION_K_H_

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#else
#include <AMReX_FArrayBox.H>
#endif

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_xslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i  ,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real dr = 2.0*(vcc(i+1,j,k,n) - vcc(i  ,j,k,n));
    amrex::Real dc = 0.5*(vcc(i+1,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_ho_xslope (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp;
    amrex::Real qm, qp, qi;
    qi = q(i,j,k,n);
    qm = q(i-1,j,k,n);
    qp = q(i+1,j,k,n);

    dlft = qm - q(i-2,j,k,n);
    drgt = qi - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qp - qi;
    drgt = q(i+2,j,k,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qi - qm;
    drgt = qp - qi;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);

    return dsgn*amrex::min(dlim, amrex::Math::abs(dtemp));

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_xslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i  ,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real dr = 2.0*(vcc(i+1,j,k,n) - vcc(i  ,j,k,n));
    amrex::Real dc = 0.5*(vcc(i+1,j,k,n) - vcc(i-1,j,k,n));
    if (edlo and i == domlo) {
        dc = (vcc(i+1,j,k,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i-1,j,k,n))/3.0;
    } else if (edhi and i == domhi) {
        dc = (4.0*vcc(i+1,j,k,n)-3.0*vcc(i,j,k,n)-vcc(i-1,j,k,n))/3.0;
    }
    amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_ho_xslope_extdir (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& q,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, dlimsh, dsgnsh;
    amrex::Real qm, qp, qi;
    qi = q(i,j,k,n);
    qm = q(i-1,j,k,n);
    qp = q(i+1,j,k,n);

    dlft = qm - q(i-2,j,k,n);
    drgt = qi - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qp - qi;
    drgt = q(i+2,j,k,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qi - qm;
    drgt = qp - qi;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);

    if (edlo and i == domlo) {
       dtemp  = -16./15.*q(i-1,j,k,n) + .5*q(i,j,k,n) + 2./3.*q(i+1,j,k,n) -  0.1*q(i+2,j,k,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i-1,j,k,n));
       drgt = 2.*(q(i+1,j,k,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgn = amrex::Math::copysign(1.e0, dtemp);
    } else if (edlo and i == domlo+1) {
       dfm  = -16./15.*q(domlo-1,j,k,n) + .5*q(domlo,j,k,n) + 2./3.*q(domlo+1,j,k,n) -  0.1*q(domlo+2,j,k,n);
       dlft = 2.*(q(domlo  ,j,k,n)-q(domlo-1,j,k,n));
       drgt = 2.*(q(domlo+1,j,k,n)-q(domlo  ,j,k,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgnsh = amrex::Math::copysign(1.e0, dfm);
       dfm = dsgnsh*amrex::min(dlimsh, amrex::Math::abs(dfm));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }

    if (edhi and i == domhi) {
       dtemp  = 16./15.*q(i+1,j,k,n) - .5*q(i,j,k,n) - 2./3.*q(i-1,j,k,n) +  0.1*q(i-2,j,k,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i-1,j,k,n));
       drgt = 2.*(q(i+1,j,k,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgn = amrex::Math::copysign(1.e0, dtemp);
    } else if (edhi and i == domhi-1) {
       dfp  = 16./15.*q(domhi+1,j,k,n) - .5*q(domhi,j,k,n) - 2./3.*q(domhi-1,j,k,n) +  0.1*q(domhi-2,j,k,n);
       dlft = 2.*(q(domhi  ,j,k,n)-q(domhi-1,j,k,n));
       drgt = 2.*(q(domhi+1,j,k,n)-q(domhi  ,j,k,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgnsh = amrex::Math::copysign(1.e0, dfp);
       dfp = dsgnsh*amrex::min(dlimsh, amrex::Math::abs(dfp));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }

    return dsgn*amrex::min(dlim, amrex::Math::abs(dtemp));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_yslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j  ,k,n) - vcc(i,j-1,k,n));
    amrex::Real dr = 2.0*(vcc(i,j+1,k,n) - vcc(i,j  ,k,n));
    amrex::Real dc = 0.5*(vcc(i,j+1,k,n) - vcc(i,j-1,k,n));
    amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_ho_yslope (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp;
    amrex::Real qm, qp, qj;
    qj = q(i,j,k,n);
    qm = q(i,j-1,k,n);
    qp = q(i,j+1,k,n);

    dlft = qm - q(i,j-2,k,n);
    drgt = qj - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qp - qj;
    drgt = q(i,j+2,k,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qj - qm;
    drgt = qp - qj;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    return dsgn*amrex::min(dlim, amrex::Math::abs(dtemp));

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_yslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j  ,k,n) - vcc(i,j-1,k,n));
    amrex::Real dr = 2.0*(vcc(i,j+1,k,n) - vcc(i,j  ,k,n));
    amrex::Real dc = 0.5*(vcc(i,j+1,k,n) - vcc(i,j-1,k,n));
    if (edlo and j == domlo) {
        dc = (vcc(i,j+1,k,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i,j-1,k,n))/3.0;
    } else if (edhi and j == domhi) {
        dc = (4.0*vcc(i,j+1,k,n)-3.0*vcc(i,j,k,n)-vcc(i,j-1,k,n))/3.0;
    }
    amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_ho_yslope_extdir (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& q,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, dlimsh,dsgnsh;
    amrex::Real qm, qp, qj;
    qj = q(i,j,k,n);
    qm = q(i,j-1,k,n);
    qp = q(i,j+1,k,n);

    dlft = qm - q(i,j-2,k,n);
    drgt = qj - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qp - qj;
    drgt = q(i,j+2,k,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qj - qm;
    drgt = qp - qj;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);

    if (edlo and j == domlo) {
       dtemp  = -16./15.*q(i,j-1,k,n) + .5*q(i,j,k,n) + 2./3.*q(i,j+1,k,n) -  0.1*q(i,j+2,k,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i,j-1,k,n));
       drgt = 2.*(q(i,j+1,k,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgn = amrex::Math::copysign(1.e0, dtemp);
    } else if (edlo and j == domlo+1) {
       dfm  = -16./15.*q(i,domlo-1,k,n) + .5*q(i,domlo,k,n) + 2./3.*q(i,domlo+1,k,n) -  0.1*q(i,domlo+2,k,n);
       dlft = 2.*(q(i  ,domlo,k,n)-q(i,domlo-1,k,n));
       drgt = 2.*(q(i,domlo+1,k,n)-q(i  ,domlo,k,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgnsh = amrex::Math::copysign(1.e0, dfm);
       dfm = dsgnsh*amrex::min(dlimsh, amrex::Math::abs(dfm));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }

    if (edhi and j == domhi) {
       dtemp  = 16./15.*q(i,j+1,k,n) - .5*q(i,j,k,n) - 2./3.*q(i,j-1,k,n) +  0.1*q(i,j-2,k,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i,j-1,k,n));
       drgt = 2.*(q(i,j+1,k,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgn = amrex::Math::copysign(1.e0, dtemp);
    } else if (edhi and j == domhi-1) {
       dfp  = 16./15.*q(i,domhi+1,k,n) - .5*q(i,domhi,k,n) - 2./3.*q(i,domhi-1,k,n) +  0.1*q(i,domhi-2,k,n);
       dlft = 2.*(q(i  ,domhi,k,n)-q(i,domhi-1,k,n));
       drgt = 2.*(q(i,domhi+1,k,n)-q(i  ,domhi,k,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgnsh = amrex::Math::copysign(1.e0, dfp);
       dfp = dsgnsh*amrex::min(dlimsh, amrex::Math::abs(dfp));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }

    return dsgn*amrex::min(dlim, amrex::Math::abs(dtemp));

}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_zslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j,k  ,n) - vcc(i,j,k-1,n));
    amrex::Real dr = 2.0*(vcc(i,j,k+1,n) - vcc(i,j,k  ,n));
    amrex::Real dc = 0.5*(vcc(i,j,k+1,n) - vcc(i,j,k-1,n));
    amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_ho_zslope (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp;
    amrex::Real qm, qp, qk;
    qk = q(i,j,k,n);
    qm = q(i,j,k-1,n);
    qp = q(i,j,k+1,n);

    dlft = qm - q(i,j,k-2,n);
    drgt = qk - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qp - qk;
    drgt = q(i,j,k+2,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qk - qm;
    drgt = qp - qk;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    return dsgn*amrex::min(dlim, amrex::Math::abs(dtemp));

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_zslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j,k  ,n) - vcc(i,j,k-1,n));
    amrex::Real dr = 2.0*(vcc(i,j,k+1,n) - vcc(i,j,k  ,n));
    amrex::Real dc = 0.5*(vcc(i,j,k+1,n) - vcc(i,j,k-1,n));
    if (edlo and k == domlo) {
        dc = (vcc(i,j,k+1,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i,j,k-1,n))/3.0;
    } else if (edhi and k == domhi) {
        dc = (4.0*vcc(i,j,k+1,n)-3.0*vcc(i,j,k,n)-vcc(i,j,k-1,n))/3.0;
    }
    amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_ho_zslope_extdir (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& q,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, dlimsh, dsgnsh;
    amrex::Real qm, qp, qk;
    qk = q(i,j,k,n);
    qm = q(i,j,k-1,n);
    qp = q(i,j,k+1,n);

    dlft = qm - q(i,j,k-2,n);
    drgt = qk - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qp - qk;
    drgt = q(i,j,k+2,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, amrex::Math::abs(dcen));

    dlft = qk - qm;
    drgt = qp - qk;
    dcen = 0.5*(dlft+drgt);
    dsgn = amrex::Math::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);

    if (edlo and k == domlo) {
       dtemp  = -16./15.*q(i,j,k-1,n) + .5*q(i,j,k,n) + 2./3.*q(i,j,k+1,n) -  0.1*q(i,j,k+2,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i,j,k-1,n));
       drgt = 2.*(q(i,j,k+1,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgn = amrex::Math::copysign(1.e0, dtemp);
    } else if (edlo and k == domlo+1) {
       dfm  = -16./15.*q(i,j,domlo-1,n) + .5*q(i,j,domlo,n) + 2./3.*q(i,j,domlo+1,n) -  0.1*q(i,j,domlo+2,n);
       dlft = 2.*(q(i  ,j,domlo,n)-q(i,j,domlo-1,n));
       drgt = 2.*(q(i,j,domlo+1,n)-q(i  ,j,domlo,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgnsh = amrex::Math::copysign(1.e0, dfm);
       dfm = dsgnsh*amrex::min(dlimsh, amrex::Math::abs(dfm));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }

    if (edhi and k == domhi) {
       dtemp  = 16./15.*q(i,j,k+1,n) - .5*q(i,j,k,n) - 2./3.*q(i,j,k-1,n) +  0.1*q(i,j,k-2,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i,j,k-1,n));
       drgt = 2.*(q(i,j,k+1,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgn = amrex::Math::copysign(1.e0, dtemp);
    } else if (edhi and k == domhi-1) {
       dfp  = 16./15.*q(i,j,domhi+1,n) - .5*q(i,j,domhi,n) - 2./3.*q(i,j,domhi-1,n) +  0.1*q(i,j,domhi-2,n);
       dlft = 2.*(q(i  ,j,domhi,n)-q(i,j,domhi-1,n));
       drgt = 2.*(q(i,j,domhi+1,n)-q(i  ,j,domhi,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt)) : 0.0;
       dsgnsh = amrex::Math::copysign(1.e0, dfp);
       dfp = dsgnsh*amrex::min(dlimsh, amrex::Math::abs(dfp));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }
    return dsgn*amrex::min(dlim, amrex::Math::abs(dtemp));


}

#ifdef AMREX_USE_EB

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
incflo_slopes_eb (int i, int j, int k, int n,
                  amrex::Array4<amrex::Real const> const& state,
                  amrex::Array4<amrex::Real const> const& ccent,
                  amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{

#if (AMREX_SPACEDIM == 3)
    constexpr int m_size = 27;
#else
    constexpr int m_size = 9;
#endif

    amrex::Real A[m_size][AMREX_SPACEDIM];
    amrex::Real du[m_size];

    {
      int lc=0;
#if (AMREX_SPACEDIM == 3)
      for(int kk(-1); kk<=1; kk++){
#else
      int kk = 0;
#endif
        for(int jj(-1); jj<=1; jj++){
          for(int ii(-1); ii<=1; ii++){

            if( flag(i,j,k).isConnected(ii,jj,kk) and
                not (ii==0 and jj==0 and kk==0)) {

            // Not multplying by dx to be consistent with how the
            // slope is stored. Also not including the global shift
            // wrt plo or i,j,k. We only need relative distance.

              A[lc][0] = ii + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
              A[lc][1] = jj + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);
#if (AMREX_SPACEDIM == 3)
              A[lc][2] = kk + ccent(i+ii,j+jj,k+kk,2) - ccent(i,j,k,2);
#endif
              du[lc] = state(i+ii,j+jj,k+kk,n) - state(i,j,k,n);

            } else {

              AMREX_D_TERM(A[lc][0] = 0.0;,
                           A[lc][1] = 0.0;,
                           A[lc][2] = 0.0;);

              du[lc] = 0.0;
            }

            lc++;
          }
        }
#if (AMREX_SPACEDIM == 3)
      }
#endif
    }

    amrex::Real AtA[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Atb[AMREX_SPACEDIM];

    for(int jj(0); jj<AMREX_SPACEDIM; ++jj){
      for(int ii(0); ii<AMREX_SPACEDIM; ++ii){
        AtA[ii][jj] = 0.0;
      }
      Atb[jj] = 0.0;
    }

    {

      for(int lc(0); lc<m_size; ++lc){
        AtA[0][0] += A[lc][0]* A[lc][0];
        AtA[0][1] += A[lc][0]* A[lc][1];
#if (AMREX_SPACEDIM == 3)
        AtA[0][2] += A[lc][0]* A[lc][2];
#endif
        AtA[1][1] += A[lc][1]* A[lc][1];
#if (AMREX_SPACEDIM == 3)
        AtA[1][2] += A[lc][1]* A[lc][2];
        AtA[2][2] += A[lc][2]* A[lc][2];
#endif

        Atb[0] += A[lc][0]*du[lc];
        Atb[1] += A[lc][1]*du[lc];
#if (AMREX_SPACEDIM == 3)
        Atb[2] += A[lc][2]*du[lc];
#endif
      }
    }

    // Fill in symmetric
    AtA[1][0] = AtA[0][1];
#if (AMREX_SPACEDIM == 3)
    AtA[2][0] = AtA[0][2];
    AtA[2][1] = AtA[1][2];
#endif

#if (AMREX_SPACEDIM == 3)
    amrex::Real detAtA =
      AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[2][1]) -
      AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
      AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);
#else
    amrex::Real detAtA =
      (AtA[0][0]*AtA[1][1])-
      (AtA[0][1]*AtA[1][0]);
#endif

    AMREX_D_TERM(amrex::Real xs = 0.0;,
                 amrex::Real ys = 0.0;,
                 amrex::Real zs = 0.0;);

    // X direction
    if(flag(i  ,j,k).isSingleValued() or
      (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
      (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0))) 
    {

#if (AMREX_SPACEDIM == 3)
      amrex::Real detAtA_x =
        Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
        AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
        AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );
#else
      amrex::Real detAtA_x =
        (Atb[0]   *AtA[1][1]) -
        (AtA[0][1]*Atb[1]); 
#endif

      // Slope at centroid of (i,j,k)
      xs = detAtA_x / detAtA;

    } else {

      amrex::Real du_xl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
      amrex::Real du_xr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
      amrex::Real du_xc = 0.5*(state(i+1,j,k,n) - state(i-1,j,k,n));

      amrex::Real xslope = amrex::min(amrex::Math::abs(du_xl), amrex::Math::abs(du_xc), amrex::Math::abs(du_xr));
      xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
      xs = (du_xc > 0.0) ? xslope : -xslope;

    }

    // Y direction
    if(flag(i,j  ,k).isSingleValued() or
      (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
      (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0))) 
    {

#if (AMREX_SPACEDIM == 3)
      amrex::Real detAtA_y =
        AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
        Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
        AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);
#else
      amrex::Real detAtA_y =
        (AtA[0][0]*Atb[1]) -
        (Atb[0] * AtA[1][0]);
#endif

      // Slope at centroid of (i,j,k)
      ys = detAtA_y / detAtA;

    } else {

      amrex::Real du_yl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
      amrex::Real du_yr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
      amrex::Real du_yc = 0.5*(state(i,j+1,k,n) - state(i,j-1,k,n));

      amrex::Real yslope = amrex::min(amrex::Math::abs(du_yl), amrex::Math::abs(du_yc), amrex::Math::abs(du_yr));
      yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
      ys = (du_yc > 0.0) ? yslope : -yslope;
    }
#if (AMREX_SPACEDIM == 3)
    // Z direction
    if(flag(i,j,k  ).isSingleValued() or
      (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
      (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1))) 
    {

      amrex::Real detAtA_z =
        AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
        AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
        Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

      zs = detAtA_z / detAtA;

    } else {

      amrex::Real du_zl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
      amrex::Real du_zr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
      amrex::Real du_zc = 0.5*(state(i,j,k+1,n) - state(i,j,k-1,n));

      amrex::Real zslope = amrex::min(amrex::Math::abs(du_zl), amrex::Math::abs(du_zc), amrex::Math::abs(du_zr));
      zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
      zs = (du_zc > 0.0) ? zslope : -zslope;
    }
#endif

   return {AMREX_D_DECL(xs,ys,zs)};
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
incflo_slopes_extdir_eb (int i, int j, int k, int n,
                         amrex::Array4<amrex::Real const> const& state,
                         amrex::Array4<amrex::Real const> const& ccent,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         AMREX_D_DECL(bool edlo_x, bool edlo_y, bool edlo_z),
                         AMREX_D_DECL(bool edhi_x, bool edhi_y, bool edhi_z), 
                         AMREX_D_DECL(int domlo_x, int domlo_y, int domlo_z),
                         AMREX_D_DECL(int domhi_x, int domhi_y, int domhi_z)) noexcept
{
    // First get EB-aware slope that doesn't know about extdir
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> eb_slope = incflo_slopes_eb (i,j,k,n,state,ccent,flag);
    AMREX_D_TERM(amrex::Real xslope = eb_slope[0];,
                 amrex::Real yslope = eb_slope[1];,
                 amrex::Real zslope = eb_slope[2];);

    // Now correct only those cells which are not affected by EB but are by extdir
    if        (edlo_x and i == domlo_x) {
        if( !( flag(i  ,j,k).isSingleValued() or
              (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
              (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0))) ) 
        {
           amrex::Real dl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
           amrex::Real dr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
           amrex::Real dc = 
              (state(i+1,j,k,n)+3.0*state(i,j,k,n)-4.0*state(i-1,j,k,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           xslope = (dc > 0.0) ? slope : -slope;
        }

    } else if (edhi_x and i == domhi_x) {
        if( !( flag(i  ,j,k).isSingleValued() or
              (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
              (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0))) )
        {
           amrex::Real dl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
           amrex::Real dr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
           amrex::Real dc = 
              (4.0*state(i+1,j,k,n)-3.0*state(i,j,k,n)-state(i-1,j,k,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           xslope = (dc > 0.0) ? slope : -slope;
        }
    }

    if        (edlo_y and j == domlo_y) {
        if( !( flag(i,j  ,k).isSingleValued() or
              (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
              (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0))) )
        {
           amrex::Real dl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
           amrex::Real dr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
           amrex::Real dc = 
              (state(i,j+1,k,n)+3.0*state(i,j,k,n)-4.0*state(i,j-1,k,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           yslope = (dc > 0.0) ? slope : -slope;
        }
    } else if (edhi_y and j == domhi_y) {
        if( !(flag(i,j  ,k).isSingleValued() or
           (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
           (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0))) )
        {
           amrex::Real dl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
           amrex::Real dr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
           amrex::Real dc = 
              (4.0*state(i,j+1,k,n)-3.0*state(i,j,k,n)-state(i,j-1,k,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           yslope = (dc > 0.0) ? slope : -slope;
        }
    }
#if (AMREX_SPACEDIM == 3)
    if        (edlo_z and k == domlo_z) {
        if( !( flag(i,j,k  ).isSingleValued() or
              (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
              (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1))) )
        {
           amrex::Real dl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
           amrex::Real dr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
           amrex::Real dc = 
              (state(i,j,k+1,n)+3.0*state(i,j,k,n)-4.0*state(i,j,k-1,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           zslope = (dc > 0.0) ? slope : -slope;
        }
    } else if (edhi_z and k == domhi_z) {
        if( !( flag(i,j,k  ).isSingleValued() or
             (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
             (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1))) )
        {
           amrex::Real dl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
           amrex::Real dr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
           amrex::Real dc = 
              (4.0*state(i,j,k+1,n)-3.0*state(i,j,k,n)-state(i,j,k-1,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           zslope = (dc > 0.0) ? slope : -slope;
        }
    }
#endif

    return {AMREX_D_DECL(xslope,yslope,zslope)};
}

#endif

}
#endif
