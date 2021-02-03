#ifdef AMREX_USE_EB
#include <Redistribution.H>

using namespace amrex;

void redistribution::flux_redistribute_eb (Box const& bx, int ncomp,
                                           Array4<Real> const& dUdt,
                                           Array4<Real const> const& dUdt_in,
                                           Array4<Real> const& scratch,
                                           Array4<EBCellFlag const> const& flag,
                                           Array4<Real const> const& vfrac,
                                           Geometry& lev_geom)
{
    const Box dbox = lev_geom.growPeriodicDomain(2);

    Array4<Real> tmp(scratch, 0);
    Array4<Real> delm(scratch, ncomp);
    Array4<Real> wgt(scratch, 2*ncomp);

    Box const& bxg1 = amrex::grow(bx,1);
    Box const& bxg2 = amrex::grow(bx,2);

    // xxxxx TODO: more weight options
    amrex::ParallelFor(bxg2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        wgt(i,j,k) = (dbox.contains(IntVect(AMREX_D_DECL(i,j,k)))) ? 1.0 : 0.0;
    });

    amrex::ParallelFor(bxg1, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isSingleValued()) {
#if 1
            Real vtot = 0.0;
            Real divnc = 0.0;
            for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
                 if ((ii != 0 or jj != 0 or kk != 0) and
                     flag(i,j,k).isConnected(ii,jj,kk) and
                    dbox.contains(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk))))
                {
                    Real vf = vfrac(i+ii,j+jj,k+kk);
                    vtot += vf;
                    divnc += vf * dUdt_in(i+ii,j+jj,k+kk,n);
                }
            }}}
            divnc /= (vtot + 1.e-80);
#else
            Real divnc = vfrac(i,j,k)*dUdt_in(i,j,k,n);
#endif
            Real optmp = (1.0-vfrac(i,j,k))*(divnc-dUdt_in(i,j,k,n));
            tmp(i,j,k,n) = optmp;
            delm(i,j,k,n) = -vfrac(i,j,k)*optmp;
        } else {
            tmp(i,j,k,n) = 0.0;
        }
    });

    amrex::ParallelFor(bxg1 & dbox, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isSingleValued()) {
            Real wtot = 0.0;
            for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
                if ((ii != 0 or jj != 0 or kk != 0) and
                    flag(i,j,k).isConnected(ii,jj,kk))
                {
                    wtot += vfrac(i+ii,j+jj,k+kk) * wgt(i+ii,j+jj,k+kk);
                }
            }}}
            wtot = 1.0/(wtot+1.e-80);

            Real dtmp = delm(i,j,k,n) * wtot;
            for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
                if ((ii != 0 or jj != 0 or kk != 0) and
                    bx.contains(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk))) and
                    flag(i,j,k).isConnected(ii,jj,kk))
                {
                    Gpu::Atomic::Add(&tmp(i+ii,j+jj,k+kk,n), dtmp*wgt(i+ii,j+jj,k+kk));
                }
            }}}
        }
    });

    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        dUdt(i,j,k,n) = dUdt_in(i,j,k,n) + tmp(i,j,k,n);
    });
}
#endif
