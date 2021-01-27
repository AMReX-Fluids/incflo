#include <MOL.H>

using namespace amrex;

void 
mol::compute_convective_rate (Box const& bx, int ncomp,
                              Array4<Real> const& dUdt,
                              AMREX_D_DECL(Array4<Real const> const& fx,
                                           Array4<Real const> const& fy,
                                           Array4<Real const> const& fz),
                              Geometry& geom)
{
    const auto dxinv = geom.InvCellSizeArray();
    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
#if (AMREX_SPACEDIM == 3)
        dUdt(i,j,k,n) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
            +           dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n))
            +           dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n));
#else
        dUdt(i,j,k,n) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
            +           dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n));
#endif
    });
}

#ifdef AMREX_USE_EB
void 
mol::compute_convective_rate_eb (Box const& bx, int ncomp,
                                 Array4<Real> const& dUdt,
                                 AMREX_D_DECL(Array4<Real const> const& fx,
                                              Array4<Real const> const& fy,
                                              Array4<Real const> const& fz),
                                 Array4<EBCellFlag const> const& flag,
                                 Array4<Real const> const& vfrac,
                                 AMREX_D_DECL(Array4<Real const> const& apx,
                                              Array4<Real const> const& apy,
                                              Array4<Real const> const& apz),
                                 Geometry& geom)
{
    const auto dxinv = geom.InvCellSizeArray();
    const Box dbox   = geom.growPeriodicDomain(2);
    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
#if (AMREX_SPACEDIM == 3)
        if (!dbox.contains(IntVect(AMREX_D_DECL(i,j,k))) or flag(i,j,k).isCovered()) {
            dUdt(i,j,k,n) = 0.0;
        } else if (flag(i,j,k).isRegular()) {
            dUdt(i,j,k,n) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
                +           dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n))
                +           dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n));
        } else {
            dUdt(i,j,k,n) = (1.0/vfrac(i,j,k)) *
                ( dxinv[0] * (apx(i,j,k)*fx(i,j,k,n) - apx(i+1,j,k)*fx(i+1,j,k,n))
                + dxinv[1] * (apy(i,j,k)*fy(i,j,k,n) - apy(i,j+1,k)*fy(i,j+1,k,n))
                + dxinv[2] * (apz(i,j,k)*fz(i,j,k,n) - apz(i,j,k+1)*fz(i,j,k+1,n)) );
        }
#else
        if (!dbox.contains(IntVect(AMREX_D_DECL(i,j,k))) or flag(i,j,k).isCovered()) {
            dUdt(i,j,k,n) = 0.0;
        } else if (flag(i,j,k).isRegular()) {
            dUdt(i,j,k,n) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
                +           dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n));
        } else {
            dUdt(i,j,k,n) = (1.0/vfrac(i,j,k)) *
                ( dxinv[0] * (apx(i,j,k)*fx(i,j,k,n) - apx(i+1,j,k)*fx(i+1,j,k,n))
                + dxinv[1] * (apy(i,j,k)*fy(i,j,k,n) - apy(i,j+1,k)*fy(i,j+1,k,n)) );
        }
#endif
    });
}
#endif
