#ifndef DERIVE_K_H_
#define DERIVE_K_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate (int i, int j, int k,
                               AMREX_D_DECL(amrex::Real idx, amrex::Real idy, amrex::Real idz),
                               amrex::Array4<amrex::Real const> const& vel) noexcept
{
    using namespace amrex;

#if (AMREX_SPACEDIM == 3)
    Real ux = amrex::Real(0.5) * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
    Real vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
    Real wx = amrex::Real(0.5) * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;

    Real uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
    Real vy = amrex::Real(0.5) * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
    Real wy = amrex::Real(0.5) * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;

    Real uz = amrex::Real(0.5) * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
    Real vz = amrex::Real(0.5) * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;
    Real wz = amrex::Real(0.5) * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz;

    return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy + amrex::Real(2.0) * wz*wz
                     + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
#else
    Real ux = amrex::Real(0.5) * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
    Real vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;

    Real uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
    Real vy = amrex::Real(0.5) * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;

    return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy
                     + (uy+vx)*(uy+vx));
#endif
}

#ifdef AMREX_USE_EB
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate_eb (int i, int j, int k,
                                  AMREX_D_DECL(amrex::Real idx, amrex::Real idy, amrex::Real idz),
                                  amrex::Array4<amrex::Real const> const& vel,
                                  amrex::EBCellFlag flag) noexcept
{
    using namespace amrex;

    if (flag.isRegular()) {
        return incflo_strainrate(i,j,k,AMREX_D_DECL(idx,idy,idz),vel);
    } else if (flag.isSingleValued()) {
        constexpr Real c0 = -1.5;
        constexpr Real c1 =  amrex::Real(2.0);
        constexpr Real c2 = -amrex::Real(0.5);
        Real ux, vx, uy, vy;
#if (AMREX_SPACEDIM == 3)
        Real wx, wy, uz, vz, wz;
#endif
        if (!flag.isConnected( 1,0,0))
        {
            // Covered cell to the right, go fish left
            ux = - (c0 * vel(i  ,j,k,0)
                  + c1 * vel(i-1,j,k,0)
                  + c2 * vel(i-2,j,k,0)) * idx;
            vx = - (c0 * vel(i  ,j,k,1)
                  + c1 * vel(i-1,j,k,1)
                  + c2 * vel(i-2,j,k,1)) * idx;
#if (AMREX_SPACEDIM == 3)
            wx = - (c0 * vel(i  ,j,k,2)
                  + c1 * vel(i-1,j,k,2)
                  + c2 * vel(i-2,j,k,2)) * idx;
#endif
        }
        else if (!flag.isConnected(-1,0,0))
        {
            // Covered cell to the left, go fish right
            ux = (c0 * vel(i  ,j,k,0)
                + c1 * vel(i+1,j,k,0)
                + c2 * vel(i+2,j,k,0)) * idx;
            vx = (c0 * vel(i  ,j,k,1)
                + c1 * vel(i+1,j,k,1)
                + c2 * vel(i+2,j,k,1)) * idx;
#if (AMREX_SPACEDIM == 3)
            wx = (c0 * vel(i  ,j,k,2)
                + c1 * vel(i+1,j,k,2)
                + c2 * vel(i+2,j,k,2)) * idx;
#endif
        }
        else
        {
            // No covered cells right or left, use standard stencil
            ux = amrex::Real(0.5) * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
            vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
#if (AMREX_SPACEDIM == 3)
            wx = amrex::Real(0.5) * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;
#endif
        }

        // Do the same in y-direction
        if (!flag.isConnected(0, 1,0))
        {
            uy = - (c0 * vel(i,j  ,k,0)
                  + c1 * vel(i,j-1,k,0)
                  + c2 * vel(i,j-2,k,0)) * idy;
            vy = - (c0 * vel(i,j  ,k,1)
                  + c1 * vel(i,j-1,k,1)
                  + c2 * vel(i,j-2,k,1)) * idy;
#if (AMREX_SPACEDIM == 3)
            wy = - (c0 * vel(i,j  ,k,2)
                  + c1 * vel(i,j-1,k,2)
                  + c2 * vel(i,j-2,k,2)) * idy;
#endif
        }
        else if (!flag.isConnected(0,-1,0))
        {
            uy = (c0 * vel(i,j  ,k,0)
                + c1 * vel(i,j+1,k,0)
                + c2 * vel(i,j+2,k,0)) * idy;
            vy = (c0 * vel(i,j  ,k,1)
                + c1 * vel(i,j+1,k,1)
                + c2 * vel(i,j+2,k,1)) * idy;
#if (AMREX_SPACEDIM == 3)
            wy = (c0 * vel(i,j  ,k,2)
                + c1 * vel(i,j+1,k,2)
                + c2 * vel(i,j+2,k,2)) * idy;
#endif
        }
        else
        {
            uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
            vy = amrex::Real(0.5) * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
#if (AMREX_SPACEDIM == 3)
            wy = amrex::Real(0.5) * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;
#endif
        }
#if (AMREX_SPACEDIM == 3)
        // Do the same in z-direction
        if (!flag.isConnected(0,0, 1))
        {
            uz = - (c0 * vel(i,j,k  ,0)
                  + c1 * vel(i,j,k-1,0)
                  + c2 * vel(i,j,k-2,0)) * idz;
            vz = - (c0 * vel(i,j,k  ,1)
                  + c1 * vel(i,j,k-1,1)
                  + c2 * vel(i,j,k-2,1)) * idz;
            wz = - (c0 * vel(i,j,k  ,2)
                  + c1 * vel(i,j,k-1,2)
                  + c2 * vel(i,j,k-2,2)) * idz;
        }
        else if (!flag.isConnected(0,0,-1))
        {
            uz = (c0 * vel(i,j,k  ,0)
                + c1 * vel(i,j,k+1,0)
                + c2 * vel(i,j,k+2,0)) * idz;
            vz = (c0 * vel(i,j,k  ,1)
                + c1 * vel(i,j,k+1,1)
                + c2 * vel(i,j,k+2,1)) * idz;
            wz = (c0 * vel(i,j,k  ,2)
                + c1 * vel(i,j,k+1,2)
                + c2 * vel(i,j,k+2,2)) * idz;
        }
        else
        {
            uz = amrex::Real(0.5) * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
            vz = amrex::Real(0.5) * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;
            wz = amrex::Real(0.5) * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz;
        }

        return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy + amrex::Real(2.0) * wz*wz
                         + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
#else
        return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy
                         + (uy+vx)*(uy+vx));
#endif
    } else {
        return 0.0;
    }
}
#endif

#endif
