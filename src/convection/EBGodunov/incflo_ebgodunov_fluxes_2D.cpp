#include <Godunov.H>
#include <EBGodunov.H>
#include <incflo_godunov_trans_bc.H>

#include <AMReX_MultiCutFab.H>
#include <AMReX_EBMultiFabUtil_2D_C.H>

using namespace amrex;

void
ebgodunov::compute_godunov_fluxes (Box const& bx, int flux_comp, int ncomp,
                                   Array4<Real      > const& fx,
                                   Array4<Real      > const& fy,
                                   Array4<Real const> const& q,
                                   Array4<Real const> const& umac,
                                   Array4<Real const> const& vmac,
                                   Array4<Real const> const& fq,
                                   Array4<Real const> const& /*divu*/,
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
    Box const& xbx = amrex::surroundingNodes(bx,0);
    Box const& ybx = amrex::surroundingNodes(bx,1);
    Box const& bxg1 = amrex::grow(bx,1);

    // Start with above and grow 1 tangentially
    Box xebx = Box(xbx).grow(1,1);
    Box yebx = Box(ybx).grow(0,1);

    const Real dx = geom.CellSize(0);
    const Real dy = geom.CellSize(1);
    Real dtdx = l_dt/dx;
    Real dtdy = l_dt/dy;

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    const auto dxinv = geom.InvCellSizeArray();

    Array4<Real> Imx = makeArray4(p, bxg1, ncomp);
    p +=         Imx.size();
    Array4<Real> Ipx = makeArray4(p, bxg1, ncomp);
    p +=         Ipx.size();
    Array4<Real> Imy = makeArray4(p, bxg1, ncomp);
    p +=         Imy.size();
    Array4<Real> Ipy = makeArray4(p, bxg1, ncomp);
    p +=         Ipy.size();

    Array4<Real> xlo = makeArray4(p, xebx, ncomp);
    p +=         xlo.size();
    Array4<Real> xhi = makeArray4(p, xebx, ncomp);
    p +=         xhi.size();
    Array4<Real> ylo = makeArray4(p, yebx, ncomp);
    p +=         ylo.size();
    Array4<Real> yhi = makeArray4(p, yebx, ncomp);
    p +=         yhi.size();

    Array4<Real> xyzlo = makeArray4(p, bxg1, ncomp);
    p +=         xyzlo.size();
    Array4<Real> xyzhi = makeArray4(p, bxg1, ncomp);
    p +=         xyzhi.size();

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

    for (int n = 0; n < ncomp; n++) 
       if (!iconserv[n]) amrex::Abort("Trying to update in non-conservative in ebgodunov");

    ebgodunov::plm_fpu_x (xebx, ncomp, Imx, Ipx, q, umac,
                          flag_arr, vfrac_arr,
                          AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                          geom, l_dt, h_bcrec, pbc, is_velocity);
    ebgodunov::plm_fpu_y (yebx, ncomp, Imy, Ipy, q, vmac, 
                          flag_arr, vfrac_arr,
                          AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                          geom, l_dt, h_bcrec, pbc, is_velocity);

    amrex::ParallelFor(
        xebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (apx(i,j,k) > 0.)
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
            } else {
                Imx(i,j,k,n) = 0.;
            }

        },
        yebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (apy(i,j,k) > 0.)
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
            } else {
                Imy(i,j,k,n) = 0.;
            }
        });

    // We can reuse the space in Ipx, Ipy and Ipz.


    //
    // Upwinding on y-faces to use as transverse terms for x-faces
    //
    Array4<Real> yzlo = makeArray4(xyzlo.dataPtr(), yebx, ncomp);
    amrex::ParallelFor(
    Box(yzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (apy(i,j,k) > 0.)
        {
            const auto bc = pbc[n];
            Real l_yzlo, l_yzhi;

            l_yzlo = ylo(i,j,k,n);
            l_yzhi = yhi(i,j,k,n);
            Real vad = vmac(i,j,k);
            Godunov_trans_ybc(i, j, k, n, q, l_yzlo, l_yzhi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, is_velocity);

            Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
            Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
            yzlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_yzhi + l_yzlo);
        } else {
            yzlo(i,j,k,n) = 0.;
        }
    });
    //
    Array4<Real> qx = makeArray4(Ipx.dataPtr(), xbx, ncomp);
    amrex::ParallelFor(xbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real stl, sth;

        if (apx(i,j,k) > 0.)
        {
            stl = xlo(i,j,k,n);
            sth = xhi(i,j,k,n);

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apy(i-1,j,k) > 0. && apy(i-1,j+1,k) > 0.)
            {
                Real quxl = (apx(i,j,k)*umac(i,j,k) - apx(i-1,j,k)*umac(i-1,j,k)) * q(i-1,j,k,n);
                stl += ( - (0.5*dtdx) * quxl
                         - (0.5*dtdy) * (apy(i-1,j+1,k)*yzlo(i-1,j+1,k  ,n)*vmac(i-1,j+1,k  )
                                        -apy(i-1,j  ,k)*yzlo(i-1,j  ,k  ,n)*vmac(i-1,j  ,k  )) ) / vfrac_arr(i-1,j,k);
                if (fq && vfrac_arr(i-1,j,k) > 0.)
                    stl += 0.5*l_dt*fq(i-1,j,k,n);
            }

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apy(i,j,k) > 0. && apy(i,j+1,k) > 0.) 
            {
                Real quxh = (apx(i+1,j,k)*umac(i+1,j,k) - apx(i,j,k)*umac(i,j,k)) * q(i,j,k,n);
                sth += ( - (0.5*dtdx) * quxh
                         - (0.5*dtdy)*(apy(i,j+1,k)*yzlo(i,j+1,k,n)*vmac(i,j+1,k)
                                      -apy(i,j  ,k)*yzlo(i,j  ,k,n)*vmac(i,j  ,k)) ) / vfrac_arr(i,j,k);
                if (fq && vfrac_arr(i  ,j,k) > 0.)
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
        fx(i,j,k,flux_comp+n) = umac(i,j,k) * qx(i,j,k,n);
    }); 

    //
    // Upwinding on x-faces to use as transverse terms for y-faces
    //
    Array4<Real> xzlo = makeArray4(xyzlo.dataPtr(), xebx, ncomp);
    amrex::ParallelFor(
    Box(xzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (apx(i,j,k) > 0.)
        {
            const auto bc = pbc[n];
            Real l_xzlo, l_xzhi;

            l_xzlo = xlo(i,j,k,n);
            l_xzhi = xhi(i,j,k,n);

            Real uad = umac(i,j,k);
            Godunov_trans_xbc(i, j, k, n, q, l_xzlo, l_xzhi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, is_velocity);

            Real st = (uad >= 0.) ? l_xzlo : l_xzhi;
            Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
            xzlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_xzhi + l_xzlo);
        } else {
            xzlo(i,j,k,n) = 0.;
        }
    });
    //

    Array4<Real> qy = makeArray4(Ipy.dataPtr(), ybx, ncomp);
    amrex::ParallelFor(ybx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real stl, sth;

        if (apy(i,j,k) > 0.)
        {
            stl = ylo(i,j,k,n);
            sth = yhi(i,j,k,n);

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apx(i,j-1,k) > 0. && apx(i+1,j-1,k) > 0.)
            {
                Real qvyl = (apy(i,j,k)*vmac(i,j,k) - apy(i,j-1,k)*vmac(i,j-1,k)) * q(i,j-1,k,n);
                stl += ( - (0.5*dtdy)*qvyl  
                         - (0.5*dtdx)*(apx(i+1,j-1,k)*xzlo(i+1,j-1,k  ,n)*umac(i+1,j-1,k  )
                                      -apx(i  ,j-1,k)*xzlo(i  ,j-1,k  ,n)*umac(i  ,j-1,k  )) ) / vfrac_arr(i,j-1,k);
                if (fq && vfrac_arr(i,j-1,k) > 0.)
                    stl += 0.5*l_dt*fq(i,j-1,k,n);
            }
    
            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apx(i,j,k) > 0. && apx(i+1,j,k) > 0.)
            {
                Real qvyh = (apy(i,j+1,k)*vmac(i,j+1,k) - apy(i,j,k)*vmac(i,j,k)) * q(i,j,k,n);
                sth += ( - (0.5*dtdy)*qvyh
                         - (0.5*dtdx)*(apx(i+1,j,k)*xzlo(i+1,j,k  ,n)*umac(i+1,j,k  )
                                      -apx(i  ,j,k)*xzlo(i  ,j,k  ,n)*umac(i  ,j,k  )) ) / vfrac_arr(i,j  ,k);
                if (fq && vfrac_arr(i,j,k) > 0.)
                    sth += 0.5*l_dt*fq(i,j,k,n);
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
        fy(i,j,k,flux_comp+n) = vmac(i,j,k) * qy(i,j,k,n);
    });
}
