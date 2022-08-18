/**
 * \file hydro_state_redistribute.cpp
 * \addtogroup Redistribution
 * @{
 *
 */

#include <hydro_redistribution.H>
#include <hydro_slope_limiter_K.H>
#if (AMREX_SPACEDIM == 2)
#include <hydro_eb_slopes_2D_K.H>
#elif (AMREX_SPACEDIM == 3)
#include <hydro_eb_slopes_3D_K.H>
#endif

using namespace amrex;


void
Redistribution::StateRedistribute ( Box const& bx, int ncomp,
                                    Array4<Real> const& U_out,
                                    Array4<Real> const& U_in,
                                    Array4<EBCellFlag const> const& flag,
                                    Array4<Real const> const& vfrac_old,
                                    Array4<Real const> const& vfrac_new,
                                    AMREX_D_DECL(Array4<Real const> const& fcx,
                                                 Array4<Real const> const& fcy,
                                                 Array4<Real const> const& fcz),
                                    Array4<Real const> const& ccent,
                                    amrex::BCRec  const* d_bcrec_ptr,
                                    Array4< int const> const& itracker,
                                    Array4<Real const> const& nrs,
                                    Array4<Real const> const& alpha,
                                    Array4<Real const> const& nbhd_vol,
                                    Array4<Real const> const& cent_hat,
                                    Geometry const& lev_geom,
                                    const int max_order)
{
    // Note that itracker has {4 in 2D, 8 in 3D} components and all are initialized to zero
    // We will add to the first component every time this cell is included in a merged neighborhood,
    //    either by merging or being merged
    //
    // In 2D, we identify the cells in the remaining three components with the following ordering
    //
    // ^  6 7 8
    // |  4   5
    // j  1 2 3
    //   i --->
    //
    // In 3D, We identify the cells in the remaining three components with the following ordering
    //
    //    at k-1   |   at k  |   at k+1
    //
    // ^  15 16 17 |  6 7 8  |  24 25 26
    // |  12 13 14 |  4   5  |  21 22 23
    // j  9  10 11 |  1 2 3  |  18 19 20
    //   i --->
    //
    // Note the first component of each of these arrays should never be used
    //
#if (AMREX_SPACEDIM == 2)
    amrex::GpuArray<int,9> imap{0,-1, 0, 1,-1, 1,-1, 0, 1};
    amrex::GpuArray<int,9> jmap{0,-1,-1,-1, 0, 0, 1, 1, 1};
    amrex::GpuArray<int,9> kmap{0, 0, 0, 0, 0, 0, 0, 0, 0};
#else
    amrex::GpuArray<int,27> imap{0,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    amrex::GpuArray<int,27> jmap{0,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    amrex::GpuArray<int,27> kmap{0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
#endif

    const Box domain = lev_geom.Domain();
    const int domain_ilo = domain.smallEnd(0);
    const int domain_ihi = domain.bigEnd(0);
    const int domain_jlo = domain.smallEnd(1);
    const int domain_jhi = domain.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);
#endif

    AMREX_D_TERM(const auto& is_periodic_x = lev_geom.isPeriodic(0);,
                 const auto& is_periodic_y = lev_geom.isPeriodic(1);,
                 const auto& is_periodic_z = lev_geom.isPeriodic(2););

    // amrex::Print() << " IN STATE_REDISTRIBUTE DOING BOX " << bx << " with ncomp " << ncomp << std::endl;
    // amrex::Print() << " Box(U_in) " << Box(U_in) << std::endl;
    // amrex::Print() << " Box(U_out) " << Box(U_out) << std::endl;

    Box const& bxg1 = amrex::grow(bx,1);
    Box const& bxg2 = amrex::grow(bx,2);
    Box const& bxg3 = amrex::grow(bx,3);

    Box domain_per_grown = domain;
    if (is_periodic_x) domain_per_grown.grow(0,2);
    if (is_periodic_y) domain_per_grown.grow(1,2);
#if (AMREX_SPACEDIM == 3)
    if (is_periodic_z) domain_per_grown.grow(2,2);
#endif

    // Solution at the centroid of my nbhd
    FArrayBox    soln_hat_fab (bxg3,ncomp);
    Array4<Real> soln_hat = soln_hat_fab.array();
    Elixir   eli_soln_hat = soln_hat_fab.elixir();

    // Show "VOLD / VNEW / UIN
    amrex::ParallelFor(bxg3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if ((vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1.) || i == 5 && j == 4){
            amrex::Print() << "VOLD / VNEW / UIN " << IntVect(i,j) << " " << 
                vfrac_old(i,j,0) << " " << vfrac_new(i,j,0) << " " << U_in(i,j,0,0);
                
            if (vfrac_old(i,j,k) > 0.)
                amrex::Print() << " " << U_in(i,j,0,0) * vfrac_new(i,j,k) / vfrac_old(i,j,k); 
               // U_in(i,j,0,0) = vfrac_new(i,j,k) / vfrac_old(i,j,k);

            amrex::Print() << std::endl;
        }
    });

    // Define Qhat (from Berger and Guliani)
    // Here we initialize soln_hat to equal U_in on all cells in bxg3 so that
    //      in the event we need to use soln_hat 3 cells out from the bx limits
    //      in a modified slope computation, we have a value of soln_hat to use.
    //      But we only modify soln_hat inside bxg2
    amrex::ParallelFor(bxg3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        for (int n = 0; n < ncomp; n++)
            soln_hat(i,j,k,n) = U_in(i,j,k,n);

        if ( ( (vfrac_new(i,j,k) > 0.0) || (vfrac_old(i,j,k) > 0.0) )
                                   && bxg2.contains(IntVect(AMREX_D_DECL(i,j,k)))
                                   && domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k)))) {

            // Start with U_in(i,j,k) itself
            for (int n = 0; n < ncomp; n++)
                soln_hat(i,j,k,n) = U_in(i,j,k,n) * alpha(i,j,k,0) * vfrac_old(i,j,k);
            
            // Add correction for newly uncovered cells
            for (int n = 0; n < ncomp; n++){
                if (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1. && vfrac_old(i,j,k) == 0.)
                    soln_hat(i,j,k,n) = vfrac_new(i,j,k); 
            }

            // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int r = i+imap[itracker(i,j,k,i_nbor)];
                int s = j+jmap[itracker(i,j,k,i_nbor)];
                int t = k+kmap[itracker(i,j,k,i_nbor)];

                if (domain_per_grown.contains(IntVect(AMREX_D_DECL(r,s,t))))
                {
                    for (int n = 0; n < ncomp; n++)
                        soln_hat(i,j,k,n) += U_in(r,s,t,n) * alpha(i,j,k,1) * vfrac_old(r,s,t) / nrs(r,s,t);
                }
            }
            for (int n = 0; n < ncomp; n++)  
                soln_hat(i,j,k,n) /= nbhd_vol(i,j,k);

            if (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1.0)
                amrex::Print() << "QHAT NBVOL " << IntVect(i,j) << " " << soln_hat(i,j,0,0) << " " <<  nbhd_vol(i,j,k) << std::endl;;
        }
    });

    amrex::ParallelFor(bxg1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac_old(i,j,k) > 0.0 || vfrac_new(i,j,k) > 0.0)
        {
            int num_nbors = itracker(i,j,k,0);

            if (itracker(i,j,k,0) == 0)
            {
                if (bx.contains(IntVect(AMREX_D_DECL(i,j,k))))
                {
                    for (int n = 0; n < ncomp; n++)
                        amrex::Gpu::Atomic::Add(&U_out(i,j,k,n),alpha(i,j,k,0)*nrs(i,j,k)*soln_hat(i,j,k,n));
                }

            } else {

                for (int n = 0; n < ncomp; n++)
                {
                    bool extdir_ilo = (d_bcrec_ptr[n].lo(0) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].lo(0) == amrex::BCType::hoextrap);
                    bool extdir_ihi = (d_bcrec_ptr[n].hi(0) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].hi(0) == amrex::BCType::hoextrap);
                    bool extdir_jlo = (d_bcrec_ptr[n].lo(1) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].lo(1) == amrex::BCType::hoextrap);
                    bool extdir_jhi = (d_bcrec_ptr[n].hi(1) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].hi(1) == amrex::BCType::hoextrap);
#if (AMREX_SPACEDIM == 3)
                    bool extdir_klo = (d_bcrec_ptr[n].lo(2) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].lo(2) == amrex::BCType::hoextrap);
                    bool extdir_khi = (d_bcrec_ptr[n].hi(2) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].hi(2) == amrex::BCType::hoextrap);
#endif
                    // Initialize so that the slope stencil goes from -1:1 in each diretion
                    int nx = 1; int ny = 1; int nz = 1;

                    // Do we have enough extent in each coordinate direction to use the 3x3x3 stencil
                    //    or do we need to enlarge it?
                    AMREX_D_TERM(Real x_max = -1.e30; Real x_min = 1.e30;,
                                 Real y_max = -1.e30; Real y_min = 1.e30;,
                                 Real z_max = -1.e30; Real z_min = 1.e30;);

                    Real slope_stencil_min_width = 0.5;
#if (AMREX_SPACEDIM == 2)
                    int kk = 0;
#elif (AMREX_SPACEDIM == 3)
                    for(int kk(-1); kk<=1; kk++)
#endif
                    {
                     for(int jj(-1); jj<=1; jj++)
                      for(int ii(-1); ii<=1; ii++)
                        if (flag(i,j,k).isConnected(ii,jj,kk))
                        {
                            int r = i+ii; int s = j+jj; int t = k+kk;

                            x_max = amrex::max(x_max, cent_hat(r,s,t,0)+static_cast<Real>(ii));
                            x_min = amrex::min(x_min, cent_hat(r,s,t,0)+static_cast<Real>(ii));
                            y_max = amrex::max(y_max, cent_hat(r,s,t,1)+static_cast<Real>(jj));
                            y_min = amrex::min(y_min, cent_hat(r,s,t,1)+static_cast<Real>(jj));
#if (AMREX_SPACEDIM == 3)
                            z_max = amrex::max(z_max, cent_hat(r,s,t,2)+static_cast<Real>(kk));
                            z_min = amrex::min(z_min, cent_hat(r,s,t,2)+static_cast<Real>(kk));
#endif
                        }
                    }
                    // If we need to grow the stencil, we let it be -nx:nx in the x-direction,
                    //    for example.   Note that nx,ny,nz are either 1 or 2
                    if ( (x_max-x_min) < slope_stencil_min_width ) nx = 2;
                    if ( (y_max-y_min) < slope_stencil_min_width ) ny = 2;
#if (AMREX_SPACEDIM == 3)
                    if ( (z_max-z_min) < slope_stencil_min_width ) nz = 2;
#endif

#if 0
                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> slopes_eb;
                    if (nx*ny*nz == 1)
                        // Compute slope using 3x3x3 stencil
                        slopes_eb = amrex_calc_slopes_extdir_eb(
                                                    i,j,k,n,soln_hat,cent_hat,vfrac_old,
                                                    AMREX_D_DECL(fcx,fcy,fcz),flag,
                                                    AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                                                    AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                                                    AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                                    AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                                    max_order);
                    else
                    {
                        // Compute slope using grown stencil (no larger than 5x5x5)
                        slopes_eb = amrex_calc_slopes_extdir_eb_grown(
                                                    i,j,k,n,AMREX_D_DECL(nx,ny,nz),
                                                    soln_hat,cent_hat,vfrac_old,
                                                    AMREX_D_DECL(fcx,fcy,fcz),flag,
                                                    AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                                                    AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                                                    AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                                    AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                                    max_order);
                    }

                    // We do the limiting separately because this limiter limits the slope based on the values
                    //    extrapolated to the cell centroid (cent_hat) locations - unlike the limiter in amrex
                    //    which bases the limiting on values extrapolated to the face centroids.
                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> lim_slope =
                        amrex_calc_centroid_limiter(i,j,k,n,soln_hat,flag,slopes_eb,cent_hat);

                    AMREX_D_TERM(lim_slope[0] *= slopes_eb[0];,
                                 lim_slope[1] *= slopes_eb[1];,
                                 lim_slope[2] *= slopes_eb[2];);
#endif

                    // Add to the cell itself
                    if (bx.contains(IntVect(AMREX_D_DECL(i,j,k))))
                    {
                        Real update = soln_hat(i,j,k,n);
                        // AMREX_D_TERM(update += lim_slope[0] * (ccent(i,j,k,0)-cent_hat(i,j,k,0));,
                        //              update += lim_slope[1] * (ccent(i,j,k,1)-cent_hat(i,j,k,1));,
                        //              update += lim_slope[2] * (ccent(i,j,k,2)-cent_hat(i,j,k,2)););
                        amrex::Gpu::Atomic::Add(&U_out(i,j,k,n),alpha(i,j,k,0)*nrs(i,j,k)*update);
                    } // if bx contains

                    // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
                    for (int i_nbor = 1; i_nbor <= num_nbors; i_nbor++)
                    {
                        int r = i+imap[itracker(i,j,k,i_nbor)];
                        int s = j+jmap[itracker(i,j,k,i_nbor)];
                        int t = k+kmap[itracker(i,j,k,i_nbor)];

                        if (bx.contains(IntVect(AMREX_D_DECL(r,s,t))))
                        {
                            Real update = soln_hat(i,j,k,n);
                            // AMREX_D_TERM(update += lim_slope[0] * (ccent(r,s,t,0)-cent_hat(i,j,k,0) + static_cast<Real>(r-i));,
                            //              update += lim_slope[1] * (ccent(r,s,t,1)-cent_hat(i,j,k,1) + static_cast<Real>(s-j));,
                            //              update += lim_slope[2] * (ccent(r,s,t,2)-cent_hat(i,j,k,2) + static_cast<Real>(t-k)););
                            amrex::Gpu::Atomic::Add(&U_out(r,s,t,n),alpha(i,j,k,1)*update);
                        } // if bx contains
                    } // i_nbor
                } // n
            } // num_nbors
        } // vfrac
    });

    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (vfrac_new(i,j,k) > 0.)
        {
            // This seems to help with a compiler issue ...
            Real denom = 1. / (nrs(i,j,k) + 1.e-40);
            U_out(i,j,k,n) *= denom;
            
            if (vfrac_old(i,j,k) < 1. && vfrac_new(i,j,k) == 1. && nrs(i,j,k) == 1){
                U_out(i,j,k,n) = U_in(i,j,k,n) * vfrac_old(i,j,k);
                amrex::Print() << "Adding new fix here: U_out" << IntVect(i,j) << ": " << U_out(i,j,k,n) << std::endl;
            }
            if (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1.) 
                amrex::Print() << "VFRAC / UOUT " << IntVect(i,j) << " " << 
                vfrac_new(i,j,k) << " " << U_out(i,j,k,n) << std::endl;
        }
        else
        {
            U_out(i,j,k,n) = 1.e40;
        }
    });

#if 0 
    //
    // This tests whether the redistribution procedure was conservative --
    //      only use if bx is the whole domain
    //
    {
      for (int n = 0; n < ncomp; n++)
      {
        Real sum1(0);
        Real sum2(0);
#if (AMREX_SPACEDIM == 2)
        int k = 0;
#else
        for (int k = bx.smallEnd(2); k <= domain.bigEnd(2); k++)
#endif
        for (int j = bx.smallEnd(1); j <= domain.bigEnd(1); j++)
        for (int i = bx.smallEnd(0); i <= domain.bigEnd(0); i++)
        {
            sum1 += vfrac_old(i,j,k)*U_in(i,j,k,n);
            sum2 += vfrac_old(i,j,k)*U_out(i,j,k,n);
        }
        if (std::abs(sum1-sum2) > 1.e-8 * sum1 && std::abs(sum1-sum2) > 1.e-8)
        {
           printf("SUMS DO NOT MATCH IN STATE REDIST : %f %f ",sum1,sum2);
           amrex::Abort();
        }
      }
    }
#endif
}
/** @} */
