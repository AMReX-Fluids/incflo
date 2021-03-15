#include <incflo_ebgodunov_corner_couple.H>
#include <incflo_godunov_trans_bc.H>
#include <Godunov.H>
#include <EBGodunov.H>

using namespace amrex;

void
ebgodunov::compute_godunov_fluxes (Box const& bx, int flux_comp, int ncomp,
                                   Array4<Real      > const& fx,
                                   Array4<Real      > const& fy,
                                   Array4<Real      > const& fz,
                                   Array4<Real const> const& q,
                                   Array4<Real const> const& umac,
                                   Array4<Real const> const& vmac,
                                   Array4<Real const> const& wmac,
                                   Array4<Real const> const& fq,
                                   Array4<Real const> const& divu,
                                   Real l_dt,
                                   Vector<BCRec> const& h_bcrec,
                                          BCRec const*  pbc,
                                   int const* iconserv,
                                      Real* p,
                                   Array4<EBCellFlag const> const& flag_arr,
                                   AMREX_D_DECL(Array4<Real const> const& apx,
                                                Array4<Real const> const& apy,
                                                Array4<Real const> const& apz),
                                   Array4<Real const> const& vfrac_arr,
                                   AMREX_D_DECL(Array4<Real const> const& fcx,
                                                Array4<Real const> const& fcy,
                                                Array4<Real const> const& fcz),
                                   Array4<Real const> const& ccent_arr,
                                   Geometry& geom,
                                   bool is_velocity )
{
    for (int n = 0; n < ncomp; n++)
       if (!iconserv[n]) amrex::Abort("Trying to update in non-conservative in ebgodunov");

    // bx is the cell-centered box on which we want to compute the advective update
    Box const& xbx = amrex::surroundingNodes(bx,0);
    Box const& ybx = amrex::surroundingNodes(bx,1);
    Box const& zbx = amrex::surroundingNodes(bx,2);

    // Start with above and grow 1 tangentially
    Box xebx = Box(xbx).grow(1,1).grow(2,1);
    Box yebx = Box(ybx).grow(0,1).grow(2,1);
    Box zebx = Box(zbx).grow(0,1).grow(1,1);

    Box xebx_g2(Box(bx).grow(1).grow(1,1).grow(2,1).surroundingNodes(0));
    Box yebx_g2(Box(bx).grow(1).grow(0,1).grow(2,1).surroundingNodes(1));
    Box zebx_g2(Box(bx).grow(1).grow(0,1).grow(1,1).surroundingNodes(2));

    Box const& bxg2 = amrex::grow(bx,2);

    const Real dx = geom.CellSize(0);
    const Real dy = geom.CellSize(1);
    const Real dz = geom.CellSize(2);
    Real dtdx = l_dt/dx;
    Real dtdy = l_dt/dy;
    Real dtdz = l_dt/dz;

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    const auto dxinv = geom.InvCellSizeArray();

    Array4<Real> Imx = makeArray4(p, bxg2, ncomp);
    p +=         Imx.size();
    Array4<Real> Ipx = makeArray4(p, bxg2, ncomp);
    p +=         Ipx.size();
    Array4<Real> Imy = makeArray4(p, bxg2, ncomp);
    p +=         Imy.size();
    Array4<Real> Ipy = makeArray4(p, bxg2, ncomp);
    p +=         Ipy.size();
    Array4<Real> Imz = makeArray4(p, bxg2, ncomp);
    p +=         Imz.size();
    Array4<Real> Ipz = makeArray4(p, bxg2, ncomp);
    p +=         Ipz.size();

    Array4<Real> xlo = makeArray4(p, xebx, ncomp);
    p +=         xlo.size();
    Array4<Real> xhi = makeArray4(p, xebx, ncomp);
    p +=         xhi.size();
    Array4<Real> ylo = makeArray4(p, yebx, ncomp);
    p +=         ylo.size();
    Array4<Real> yhi = makeArray4(p, yebx, ncomp);
    p +=         yhi.size();
    Array4<Real> zlo = makeArray4(p, zebx, ncomp);
    p +=         zlo.size();
    Array4<Real> zhi = makeArray4(p, zebx, ncomp);
    p +=         zhi.size();

    Array4<Real> xyzlo = makeArray4(p, bxg2, ncomp);
    p +=         xyzlo.size();
    Array4<Real> xyzhi = makeArray4(p, bxg2, ncomp);
    p +=         xyzhi.size();

    for (int n = 0; n < ncomp; n++)
       if (!iconserv[n]) amrex::Abort("Trying to update in non-conservative in ebgodunov");

    // Initialize this way out of an abundance of paranoia
    amrex::ParallelFor(
        Box(Imx), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           Imx(i,j,k,n) = i*1.e10 + j*1.e20 + k*1.30 + n*1.e4;     
        });
    amrex::ParallelFor(
        Box(Imy), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           Imy(i,j,k,n) = i*1.e10 + j*1.e20 + k*1.30 + n*1.e4;     
        });
    amrex::ParallelFor(
        Box(Imz), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           Imz(i,j,k,n) = i*1.e10 + j*1.e20 + k*1.30 + n*1.e4;     
        });

    ebgodunov::plm_fpu_x (xebx, ncomp, Imx, Ipx, q, umac,
                          flag_arr, vfrac_arr,
                          AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                          geom, l_dt, h_bcrec, pbc, is_velocity);
    ebgodunov::plm_fpu_y (yebx, ncomp, Imy, Ipy, q, vmac,
                          flag_arr, vfrac_arr,
                          AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                          geom, l_dt, h_bcrec, pbc, is_velocity);
    ebgodunov::plm_fpu_z (zebx, ncomp, Imz, Ipz, q, wmac,
                          flag_arr, vfrac_arr,
                          AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                          geom, l_dt, h_bcrec, pbc, is_velocity);

    amrex::ParallelFor(
        xebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipx(i-1,j,k,n);
            Real hi = Imx(i  ,j,k,n);

            Real uad = umac(i,j,k);

            auto bc = pbc[n];

            Godunov_trans_xbc(i, j, k, n, q, lo, hi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, is_velocity);

            xlo(i,j,k,n) = lo;
            xhi(i,j,k,n) = hi;

            Real st = (uad >= 0.) ? lo : hi;
            Real fux = (amrex::Math::abs(uad) < small_vel)? 0. : 1.;
            Imx(i,j,k,n) = fux*st + (1. - fux)*0.5*(hi + lo);
        },
        yebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipy(i,j-1,k,n);
            Real hi = Imy(i,j  ,k,n);

            Real vad = vmac(i,j,k);

            auto bc = pbc[n];

            Godunov_trans_ybc(i, j, k, n, q, lo, hi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, is_velocity);

            ylo(i,j,k,n) = lo;
            yhi(i,j,k,n) = hi;

            Real st = (vad >= 0.) ? lo : hi;
            Real fuy = (amrex::Math::abs(vad) < small_vel)? 0. : 1.;
            Imy(i,j,k,n) = fuy*st + (1. - fuy)*0.5*(hi + lo);
        },
        zebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipz(i,j,k-1,n);
            Real hi = Imz(i,j,k  ,n);

            Real wad = wmac(i,j,k);

            auto bc = pbc[n];

            Godunov_trans_zbc(i, j, k, n, q, lo, hi, bc.lo(2), bc.hi(2), dlo.z, dhi.z, is_velocity);

            zlo(i,j,k,n) = lo;
            zhi(i,j,k,n) = hi;

            Real st = (wad >= 0.) ? lo : hi;
            Real fuz = (amrex::Math::abs(wad) < small_vel) ? 0. : 1.;
            Imz(i,j,k,n) = fuz*st + (1. - fuz)*0.5*(hi + lo);
        });

    Array4<Real> xedge = Imx;
    Array4<Real> yedge = Imy;
    Array4<Real> zedge = Imz;

    //
    // Grow in x-direction
    //
    Box const& xbxtmp = amrex::grow(bx,0,1);
    Array4<Real> yzlo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(xbxtmp,1), ncomp);
    Array4<Real> zylo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(xbxtmp,2), ncomp);

    amrex::ParallelFor(
    Box(zylo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_zylo, l_zyhi;
        EBGodunov_corner_couple_zy(l_zylo, l_zyhi,
                                   i, j, k, n, l_dt, dy, true, 
                                   zlo(i,j,k,n), zhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, vmac, yedge);

        Real wad = wmac(i,j,k);
        Godunov_trans_zbc(i, j, k, n, q, l_zylo, l_zyhi, bc.lo(2), bc.hi(2), dlo.z, dhi.z, is_velocity);

        Real st = (wad >= 0.) ? l_zylo : l_zyhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
        zylo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_zyhi + l_zylo);
    },
    Box(yzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_yzlo, l_yzhi;
        EBGodunov_corner_couple_yz(l_yzlo, l_yzhi,
                                   i, j, k, n, l_dt, dz, true, 
                                   ylo(i,j,k,n), yhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, wmac, zedge);

        Real vad = vmac(i,j,k);
        Godunov_trans_ybc(i, j, k, n, q, l_yzlo, l_yzhi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, is_velocity);

        Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yzlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_yzhi + l_yzlo);
    });
    //

    Array4<Real> qx = makeArray4(Ipx.dataPtr(), xbx, ncomp);
    amrex::ParallelFor(xbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (apx(i,j,k) > 0.)
        {
            Real stl = xlo(i,j,k,n); 
            Real sth = xhi(i,j,k,n); 

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apy(i-1,j+1,k) > 0. && apy(i-1,j  ,k) > 0. && apz(i-1,j,k+1) > 0. && apz(i-1,j,k) > 0.)
            {
                Real quxl = (apx(i,j,k)*umac(i,j,k) - apx(i-1,j,k)*umac(i-1,j,k)) * q(i-1,j,k,n);
                stl += ( - (0.5*dtdx) * quxl
                         - (0.5*dtdy)*(apy(i-1,j+1,k  )*yzlo(i-1,j+1,k  ,n)*vmac(i-1,j+1,k  )
                                     - apy(i-1,j  ,k  )*yzlo(i-1,j  ,k  ,n)*vmac(i-1,j  ,k  ))
                         - (0.5*dtdz)*(apz(i-1,j  ,k+1)*zylo(i-1,j  ,k+1,n)*wmac(i-1,j  ,k+1)
                                     - apz(i-1,j  ,k  )*zylo(i-1,j  ,k  ,n)*wmac(i-1,j  ,k  )) ) / vfrac_arr(i-1,j,k);
                if (fq && vfrac_arr(i-1,j,k) > 0.)
                    stl += 0.5*l_dt*fq(i-1,j,k,n);
            }

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apy(i,j+1,k) > 0. && apy(i,j  ,k) > 0. && apz(i,j,k+1) > 0. && apz(i,j,k) > 0.)
            {
                Real quxh = (apx(i+1,j,k)*umac(i+1,j,k) - apx(i,j,k)*umac(i,j,k)) * q(i,j,k,n);
                sth += ( - (0.5*dtdx) * quxh
                         - (0.5*dtdy)*(apy(i,j+1,k  )*yzlo(i,j+1,k  ,n)*vmac(i,j+1,k  )
                                     - apy(i,j  ,k  )*yzlo(i,j  ,k  ,n)*vmac(i,j  ,k  ))
                         - (0.5*dtdz)*(apz(i,j  ,k+1)*zylo(i,j  ,k+1,n)*wmac(i,j  ,k+1)
                                     - apz(i,j  ,k  )*zylo(i,j  ,k  ,n)*wmac(i,j  ,k  )) ) / vfrac_arr(i,j,k);
                if (fq && vfrac_arr(i,j,k) > 0.)
                    sth += 0.5*l_dt*fq(i  ,j,k,n);
            }

            auto bc = pbc[n];
            Godunov_cc_xbc_lo(i, j, k, n, q, stl, sth, bc.lo(0), dlo.x, is_velocity);
            Godunov_cc_xbc_hi(i, j, k, n, q, stl, sth, bc.hi(0), dhi.x, is_velocity);

            Real temp = (umac(i,j,k) >= 0.) ? stl : sth;
            temp = (amrex::Math::abs(umac(i,j,k)) < small_vel) ? 0.5*(stl + sth) : temp;
            qx(i,j,k,n) = temp;

        } else {
             qx(i,j,k,n) = 0.;
        }
        fx(i,j,k,n) = umac(i,j,k) * qx(i,j,k,n);
    });

    //
    // Grow in y-direction
    //
    Box const& ybxtmp = amrex::grow(bx,1,1);
    Array4<Real> xzlo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(ybxtmp,0), ncomp);
    Array4<Real> zxlo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(ybxtmp,2), ncomp);

    amrex::ParallelFor(
    Box(xzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_xzlo, l_xzhi;
        EBGodunov_corner_couple_xz(l_xzlo, l_xzhi,
                                   i, j, k, n, l_dt, dz, true, 
                                   xlo(i,j,k,n),  xhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, wmac, zedge);

        Real uad = umac(i,j,k);
        Godunov_trans_xbc(i, j, k, n, q, l_xzlo, l_xzhi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, is_velocity);

        Real st = (uad >= 0.) ? l_xzlo : l_xzhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xzlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_xzhi + l_xzlo);
    },
    Box(zxlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_zxlo, l_zxhi;
        EBGodunov_corner_couple_zx(l_zxlo, l_zxhi,
                                   i, j, k, n, l_dt, dx, true, 
                                   zlo(i,j,k,n), zhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, umac, xedge);

        Real wad = wmac(i,j,k);
        Godunov_trans_zbc(i, j, k, n, q, l_zxlo, l_zxhi, bc.lo(2), bc.hi(2), dlo.z, dhi.z, is_velocity);

        Real st = (wad >= 0.) ? l_zxlo : l_zxhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
        zxlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_zxhi + l_zxlo);
    });
    //

    Array4<Real> qy = makeArray4(Ipy.dataPtr(), ybx, ncomp);
    amrex::ParallelFor(ybx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (apy(i,j,k) > 0.)
        {
            Real stl = ylo(i,j,k,n);
            Real sth = yhi(i,j,k,n);

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apx(i+1,j-1,k) > 0. && apx(i,j-1,k) > 0. && apz(i,j-1,k+1) > 0. && apz(i,j-1,k) > 0.)
            {
                Real quyl = (apy(i,j,k)*vmac(i,j,k) - apy(i,j-1,k)*vmac(i,j-1,k)) * q(i,j-1,k,n);
                stl += ( - (0.5*dtdy) * quyl
                         - (0.5*dtdx)*(apx(i+1,j-1,k  )*xzlo(i+1,j-1,k  ,n)*umac(i+1,j-1,k  )
                                     - apx(i  ,j-1,k  )*xzlo(i  ,j-1,k  ,n)*umac(i  ,j-1,k  ))
                         - (0.5*dtdz)*(apz(i  ,j-1,k+1)*zxlo(i  ,j-1,k+1,n)*wmac(i  ,j-1,k+1)
                                     - apz(i  ,j-1,k  )*zxlo(i  ,j-1,k  ,n)*wmac(i  ,j-1,k  )) ) / vfrac_arr(i,j-1,k);
                if (fq && vfrac_arr(i,j-1,k) > 0.)
                    stl += 0.5*l_dt*fq(i,j-1,k,n);
            }

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apx(i+1,j,k) > 0. && apx(i,j,k) > 0. && apz(i,j,k+1) > 0. && apz(i,j,k) > 0.)
            {
                Real quyh = (apy(i,j+1,k)*vmac(i,j+1,k) - apy(i,j,k)*vmac(i,j,k)) * q(i,j,k,n);
                sth += ( - (0.5*dtdy) * quyh
                         - (0.5*dtdx)*(apx(i+1,j,k  )*xzlo(i+1,j,k  ,n)*umac(i+1,j,k  )
                                     - apx(i  ,j,k  )*xzlo(i  ,j,k  ,n)*umac(i  ,j,k  ))
                         - (0.5*dtdz)*(apz(i  ,j,k+1)*zxlo(i  ,j,k+1,n)*wmac(i  ,j,k+1)
                                     - apz(i  ,j,k  )*zxlo(i  ,j,k  ,n)*wmac(i  ,j,k  )) ) / vfrac_arr(i,j,k);
                if (fq && vfrac_arr(i,j  ,k) > 0.)
                    sth += 0.5*l_dt*fq(i,j  ,k,n);
            }

            auto bc = pbc[n];
            Godunov_cc_ybc_lo(i, j, k, n, q, stl, sth, bc.lo(1), dlo.y, is_velocity);
            Godunov_cc_ybc_hi(i, j, k, n, q, stl, sth, bc.hi(1), dhi.y, is_velocity);

            Real temp = (vmac(i,j,k) >= 0.) ? stl : sth;
            temp = (amrex::Math::abs(vmac(i,j,k)) < small_vel) ? 0.5*(stl + sth) : temp;
            qy(i,j,k,n) = temp;
        } else {
            qy(i,j,k,n) = 0.;
        }
        fy(i,j,k,n) = vmac(i,j,k) * qy(i,j,k,n);
    });

    //
    // z-direction
    //
    Box const& zbxtmp = amrex::grow(bx,2,1);
    Array4<Real> xylo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(zbxtmp,0), ncomp);
    Array4<Real> yxlo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(zbxtmp,1), ncomp);
    amrex::ParallelFor(
    Box(xylo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_xylo, l_xyhi;
        EBGodunov_corner_couple_xy(l_xylo, l_xyhi,
                                   i, j, k, n, l_dt, dy, true, 
                                   xlo(i,j,k,n), xhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, vmac, yedge);

        Real uad = umac(i,j,k);
        Godunov_trans_xbc(i, j, k, n, q, l_xylo, l_xyhi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, is_velocity);

        Real st = (uad >= 0.) ? l_xylo : l_xyhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xylo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_xyhi + l_xylo);
    },
    Box(yxlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_yxlo, l_yxhi;
        EBGodunov_corner_couple_yx(l_yxlo, l_yxhi,
                                   i, j, k, n, l_dt, dx, true, 
                                   ylo(i,j,k,n), yhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, umac, xedge);

        Real vad = vmac(i,j,k);
        Godunov_trans_ybc(i, j, k, n, q, l_yxlo, l_yxhi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, is_velocity);

        Real st = (vad >= 0.) ? l_yxlo : l_yxhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yxlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_yxhi + l_yxlo);
    });
    //
    Array4<Real> qz = makeArray4(Ipz.dataPtr(), zbx, ncomp);
    amrex::ParallelFor(zbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (apz(i,j,k) > 0.)
        {
            Real stl = zlo(i,j,k,n);
            Real sth = zhi(i,j,k,n);

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apx(i+1,j,k-1) > 0. && apx(i,j,k-1) > 0. && apy(i,j+1,k-1) > 0. && apy(i,j,k-1) > 0.)
            {
                Real quzl = (apz(i,j,k)*wmac(i,j,k) - apz(i,j,k-1)*wmac(i,j,k-1)) * q(i,j,k-1,n);
                stl += ( - (0.5*dtdz) * quzl
                         - (0.5*dtdx)*(apx(i+1,j  ,k-1)*xylo(i+1,j  ,k-1,n)*umac(i+1,j  ,k-1)
                                      -apx(i  ,j  ,k-1)*xylo(i  ,j  ,k-1,n)*umac(i  ,j  ,k-1))
                         - (0.5*dtdy)*(apy(i  ,j+1,k-1)*yxlo(i  ,j+1,k-1,n)*vmac(i  ,j+1,k-1)
                                      -apy(i  ,j  ,k-1)*yxlo(i  ,j  ,k-1,n)*vmac(i  ,j  ,k-1)) ) / vfrac_arr(i,j,k-1);
                if (fq && vfrac_arr(i,j,k-1) > 0.)
                    stl += 0.5*l_dt*fq(i,j,k-1,n);
            }

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apx(i+1,j,k) > 0. && apx(i,j,k) > 0. && apy(i,j+1,k) > 0. && apy(i,j,k) > 0.)
            {
                Real quzh = (apz(i,j,k+1)*wmac(i,j,k+1) - apz(i,j,k)*wmac(i,j,k)) * q(i,j,k,n);
                sth += ( - (0.5*dtdz) * quzh
                         - (0.5*dtdx)*(apx(i+1,j  ,k)*xylo(i+1,j  ,k,n)*umac(i+1,j  ,k)
                                      -apx(i  ,j  ,k)*xylo(i  ,j  ,k,n)*umac(i  ,j  ,k))
                         - (0.5*dtdy)*(apy(i  ,j+1,k)*yxlo(i  ,j+1,k,n)*vmac(i  ,j+1,k)
                                      -apy(i  ,j  ,k)*yxlo(i  ,j  ,k,n)*vmac(i  ,j  ,k)) ) / vfrac_arr(i,j,k);

                if (fq && vfrac_arr(i,j,k) > 0.)
                    sth += 0.5*l_dt*fq(i,j,k,n);
            }

            auto bc = pbc[n];
            Godunov_cc_zbc_lo(i, j, k, n, q, stl, sth, bc.lo(2),  dlo.z, is_velocity);
            Godunov_cc_zbc_hi(i, j, k, n, q, stl, sth, bc.hi(2),  dhi.z, is_velocity);

            Real temp = (wmac(i,j,k) >= 0.) ? stl : sth;
            temp = (amrex::Math::abs(wmac(i,j,k)) < small_vel) ? 0.5*(stl + sth) : temp;
            qz(i,j,k,n) = temp;

        } else {
            qz(i,j,k,n) = 0.;
        }
        fz(i,j,k,n) = wmac(i,j,k) * qz(i,j,k,n);
    });
}
