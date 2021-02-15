#include <incflo_godunov_trans_bc.H>

#include <Godunov.H>
#include <EBGodunov.H>

using namespace amrex;

void ebgodunov::make_trans_velocities (AMREX_D_DECL(Box const& xbx, 
                                                    Box const& ybx, 
                                                    Box const& zbx),
                                       AMREX_D_DECL(Array4<Real> const& u_ad,
                                                    Array4<Real> const& v_ad,
                                                    Array4<Real> const& w_ad),
                                       AMREX_D_DECL(Array4<Real const> const& Imx,
                                                    Array4<Real const> const& Imy,
                                                    Array4<Real const> const& Imz),
                                       AMREX_D_DECL(Array4<Real const> const& Ipx,
                                                    Array4<Real const> const& Ipy,
                                                    Array4<Real const> const& Ipz),
                                       Array4<Real const> const& vel,
                                       Array4<EBCellFlag const> const& flag,
                                       const Box& domain,
                                       BCRec  const* pbc)
{
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    amrex::ParallelFor(AMREX_D_DECL(xbx, ybx, zbx),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about x-velocity on x-faces here
        if (flag(i,j,k).isConnected(-1,0,0))
        {
            constexpr int n = 0;

            Real lo = Ipx(i-1,j,k,n);
            Real hi = Imx(i  ,j,k,n);

            auto bc = pbc[n];
            Godunov_trans_xbc(i, j, k, n, vel, lo, hi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, true);

            Real st = ( (lo+hi) >= 0.) ? lo : hi;
            bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
            u_ad(i,j,k) = ltm ? 0. : st;
        } else {
            u_ad(i,j,k) = 0.;
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about y-velocity on y-faces here
        if (flag(i,j,k).isConnected(0,-1,0))
        {
            constexpr int n = 1;

            Real lo = Ipy(i,j-1,k,n);
            Real hi = Imy(i,j  ,k,n);

            auto bc = pbc[n];
            Godunov_trans_ybc(i, j, k, n, vel, lo, hi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, true);

            Real st = ( (lo+hi) >= 0.) ? lo : hi;
            bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
            v_ad(i,j,k) = ltm ? 0. : st;
        } else {
            v_ad(i,j,k) = 0.;
        }
#if (AMREX_SPACEDIM == 3)
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(0,0,-1))
        {
            // We only care about z-velocity on z-faces here
            constexpr int n = 2;

            Real lo = Ipz(i,j,k-1,n);
            Real hi = Imz(i,j,k  ,n);

            auto bc = pbc[n];
            Godunov_trans_zbc(i, j, k, n, vel, lo, hi, bc.lo(2), bc.hi(2), dlo.z, dhi.z, true);

            Real st = ( (lo+hi) >= 0.) ? lo : hi;
            bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
            w_ad(i,j,k) = ltm ? 0. : st;
        } else {
            w_ad(i,j,k) = 0.;
        } 
#endif
    });
}
