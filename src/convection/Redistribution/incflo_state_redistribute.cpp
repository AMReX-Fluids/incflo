#ifdef AMREX_USE_EB
#include <Redistribution.H>
#include <AMReX_EB_slopes_K.H>

using namespace amrex;

void 
redistribution::state_redistribute ( Box const& bx, int ncomp,
                                     Array4<Real> const& dUdt,
                                     Array4<Real> const& dUdt_in,
                                     Array4<EBCellFlag const> const& flag,
                                     AMREX_D_DECL(Array4<Real const> const& apx,
                                                  Array4<Real const> const& apy,
                                                  Array4<Real const> const& apz),
                                     Array4<Real const> const& vfrac,
                                     AMREX_D_DECL(Array4<Real const> const& fcx,
                                                  Array4<Real const> const& fcy,
                                                  Array4<Real const> const& fcz),
                                     Array4<Real const> const& ccent,
                                     Geometry& lev_geom)
{
    // We identify the cells with the following ordering in 2D
    //
    // ^  6 7 8
    // |  3 4 5
    // j  0 1 2
    //   i --->
    //
    // and in 3D:
    //
    //    at k-1   |   at k     |   at k+1 
    // 
    // ^   6  7  8 |  15 16 17  |  24 25 26
    // |   3  4  5 |  12 13 14  |  21 22 23
    // j   0  1  2 |   9 10 11  |  18 19 20
    //   i --->
    // 

    const Box domain = lev_geom.Domain();

    AMREX_D_TERM(const auto& is_periodic_x = lev_geom.isPeriodic(0);,
                 const auto& is_periodic_y = lev_geom.isPeriodic(1);,
                 const auto& is_periodic_z = lev_geom.isPeriodic(2););

//  amrex::Print() << " IN STATE_REDISTRIBUTE DOING BOX " << bx << " with ncomp " << ncomp << std::endl;

    Box const& bxg1 = amrex::grow(bx,1);
    Box const& bxg2 = amrex::grow(bx,2);

    Box bx_per_grown = bx;
    if (is_periodic_x) bx_per_grown.grow(0,1);
    if (is_periodic_y) bx_per_grown.grow(1,1);
#if (AMREX_SPACEDIM == 3)
    if (is_periodic_z) bx_per_grown.grow(2,1);
#endif

    // Set to 1 if cell is in my nbhd, otherwise 0
#if (AMREX_SPACEDIM == 2)
    IArrayBox nbor_fab(bxg1,9);
#else
    IArrayBox nbor_fab(bxg1,27);
#endif

    // How many nbhds is this cell in
    FArrayBox nrs_fab       (bxg2,1);

    // Total volume of all cells in my nbhd
    FArrayBox nbhd_vol_fab  (bxg2,1);

    // Centroid of my nbhd
    FArrayBox cent_hat_fab  (bxg2,AMREX_SPACEDIM);

    // Solution at the centroid of my nbhd
    FArrayBox soln_hat_fab  (bxg2,ncomp);

    Array4<int>  nbor     = nbor_fab.array();
    Array4<Real> nbhd_vol = nbhd_vol_fab.array();
    Array4<Real> nrs      = nrs_fab.array();
    Array4<Real> soln_hat = soln_hat_fab.array();
    Array4<Real> cent_hat = cent_hat_fab.array();

    amrex::ParallelFor(bxg1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
#if (AMREX_SPACEDIM == 2)
	for (int n = 0; n < 9; n++)
#else
	for (int n = 0; n < 27; n++)
#endif
            nbor(i,j,k,n) = 0;
    });
    amrex::ParallelFor(bxg2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        nrs(i,j,k) = 0.;
        nbhd_vol(i,j,k) = 0.;
	for (int n = 0; n < AMREX_SPACEDIM; n++)
	{
            cent_hat(i,j,k,n) = 0.;
	}
	for (int n = 0; n < ncomp; n++)
	{
            soln_hat(i,j,k,n) = 0.;
	}
    });

    // It is essential that only one of these be true;
    bool vertical_only   = false;
    bool horizontal_only = true;
    bool all_four        = false;

    amrex::ParallelFor(bxg1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (!flag(i,j,k).isCovered() && bx_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k))))
        {
          // Always include the cell itself
          nbor(i,j,k,4) = 1;

          if (vfrac(i,j,k) < 0.5)
          {
            // We only include cells into a neighborhood if they are in the interior
            //    or in periodic ghost cells
            AMREX_D_TERM(bool allow_lo_x = (i > domain.smallEnd(0) || is_periodic_x);,
                         bool allow_lo_y = (j > domain.smallEnd(1) || is_periodic_y);,
                         bool allow_lo_z = (k > domain.smallEnd(2) || is_periodic_z););
            AMREX_D_TERM(bool allow_hi_x = (i < domain.bigEnd(0)   || is_periodic_x);,
                         bool allow_hi_y = (j < domain.bigEnd(1)   || is_periodic_y);,
                         bool allow_hi_z = (k < domain.bigEnd(2)   || is_periodic_z););

            if (apx(i,j,k) > 0. && allow_lo_x)
            {
                if (all_four)
                {
                    if (fcx(i,j,k,0) <= 0. && allow_lo_y)
                    {
                        if (vfrac(i-1,j-1,k) > 0.) nbor(i,j,k,0) = 1;
                        if (vfrac(i  ,j-1,k) > 0.) nbor(i,j,k,1) = 1;
                    } else if (allow_hi_y) {
                        if (vfrac(i-1,j+1,k) > 0.) nbor(i,j,k,6) = 1;
                        if (vfrac(i  ,j+1,k) > 0.) nbor(i,j,k,7) = 1;
                    }
                }
                if (all_four || horizontal_only)
                     if (vfrac(i-1,j  ,k) > 0.) nbor(i,j,k,3) = 1;
            }

            if (apx(i+1,j,k) > 0. && allow_hi_x)
            {
                if (all_four)
                {
                    if (fcx(i+1,j,k,0) <= 0. && allow_lo_y)
                    {
                        if (vfrac(i+1,j-1,k) > 0.) nbor(i,j,k,2) = 1;
                        if (vfrac(i  ,j-1,k) > 0.) nbor(i,j,k,1) = 1;
                    } else if (allow_hi_y) {
                        if (vfrac(i+1,j+1,k) > 0.) nbor(i,j,k,8) = 1;
                        if (vfrac(i  ,j+1,k) > 0.) nbor(i,j,k,7) = 1;
                    }
                }
                if (all_four || horizontal_only)
                     if (vfrac(i+1,j  ,k) > 0.) nbor(i,j,k,5) = 1;
            }

            if (apy(i,j,k) > 0. && allow_lo_y)
            {
                if (all_four)
                {
                    if (fcy(i,j,k,0) <= 0. && allow_lo_x) {
                        if (vfrac(i-1,j  ,k) > 0.) nbor(i,j,k,3) = 1;
                        if (vfrac(i-1,j-1,k) > 0.) nbor(i,j,k,0) = 1;
                    } else if (allow_hi_x) {
                        if (vfrac(i+1,j  ,k) > 0.) nbor(i,j,k,5) = 1;
                        if (vfrac(i+1,j-1,k) > 0.) nbor(i,j,k,2) = 1;
                    }
                }
                if (all_four || vertical_only)
                     if (vfrac(i  ,j-1,k) > 0.) nbor(i,j,k,1) = 1;
            }

            if (apy(i,j+1,k) > 0. && allow_hi_y)
            {
                if (all_four)
                {
                    if (fcy(i,j+1,k,0) <= 0. && allow_lo_x)
                    {
                        if (vfrac(i-1,j  ,k) > 0.) nbor(i,j,k,3) = 1;
                        if (vfrac(i-1,j+1,k) > 0.) nbor(i,j,k,6) = 1;
                    } else if (allow_hi_x) {
                        if (vfrac(i+1,j  ,k) > 0.) nbor(i,j,k,5) = 1;
                        if (vfrac(i+1,j+1,k) > 0.) nbor(i,j,k,8) = 1;
                    }
                }
                if (all_four || vertical_only)
                     if (vfrac(i  ,j+1,k) > 0.) nbor(i,j,k,7) = 1;
            }
        } // if (vfrac < 0.5)

#if (AMREX_SPACEDIM == 2)
        int kk = 0;
        for (int jj = -1; jj <= 1; jj++)  
        for (int ii = -1; ii <= 1; ii++)  
        {
            int index = (jj+1)*3 + (ii+1);
#else
        for (int kk = -1; kk <= 1; kk++)  
        for (int jj = -1; jj <= 1; jj++)  
        for (int ii = -1; ii <= 1; ii++)  
        {
            int index = (kk+1)*9 + (jj+1)*3 + (ii+1);
#endif
            if (nbor(i,j,k,index) == 1)
            {
                int r = i+ii; int s = j+jj; int t = k+kk;
                nrs(r,s,t) += 1.;
            }
        }

#if 0
        for (int index = 0; index < 9; index++)
        {
            if (index != 4 and i > 120)
                if (nbor(i,j,k,index) == 1) amrex::Print() << IntVect(i,j) << " HAS NBOR IN DIR " << index << " and vfrac " <<
                     vfrac(i,j,k) << std::endl;
        }
#endif
      }
    });

    amrex::ParallelFor(bxg1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (!flag(i,j,k).isCovered())
        {
            Real unwted_vol = 0.; // This is used as a diagnostic to make sure we don't miss any small cells
            if ( ( (i >= domain.smallEnd(0) && i <= domain.bigEnd(0)) || is_periodic_x ) &&
#if (AMREX_SPACEDIM == 3)
                 ( (k >= domain.smallEnd(2) && k <= domain.bigEnd(2)) || is_periodic_z ) &&
#endif
                 ( (j >= domain.smallEnd(1) && j <= domain.bigEnd(1)) || is_periodic_y ) )
            {
#if (AMREX_SPACEDIM == 2)
                int kk = 0;
                for (int jj = -1; jj <= 1; jj++)  
                for (int ii = -1; ii <= 1; ii++)  
                {
                    int index = (jj+1)*3 + (ii+1);
#else
                for (int kk = -1; kk <= 1; kk++)  
                for (int jj = -1; jj <= 1; jj++)  
                for (int ii = -1; ii <= 1; ii++)  
                {
                    int index = (kk+1)*9 + (jj+1)*3 + (ii+1);
#endif
                    if (nbor(i,j,k,index) == 1)
                    {
                        int r = i+ii; int s = j+jj; int t = k+kk;
                        nbhd_vol(i,j,k) += vfrac(r,s,t) / nrs(r,s,t);
                        unwted_vol += vfrac(r,s,t);
                    }
                }
            }
            if (unwted_vol < 0.5 && bx_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k))))
            {
                amrex::Print() << "NBHD VOL STILL TOO LOW " << IntVect(AMREX_D_DECL(i,j,k)) << " " << nbhd_vol(i,j,k) << std::endl;
            }
        }
    });

    // Define xhat,yhat,zhat (from Berger and Guliani) 
    amrex::ParallelFor(bxg2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac(i,j,k) > 0.5)
        {
            AMREX_D_TERM(cent_hat(i,j,k,0) = ccent(i,j,k,0);,
                         cent_hat(i,j,k,1) = ccent(i,j,k,1);,
                         cent_hat(i,j,k,2) = ccent(i,j,k,2););

        } else if (vfrac(i,j,k) > 0.0 &&
           ((i >= domain.smallEnd(0) && i <= domain.bigEnd(0)) || is_periodic_x ) &&
#if (AMREX_SPACEDIM == 3)
           ((k >= domain.smallEnd(2) && k <= domain.bigEnd(2)) || is_periodic_z ) &&
#endif
           ((j >= domain.smallEnd(1) && j <= domain.bigEnd(1)) || is_periodic_y ) )
        {
#if (AMREX_SPACEDIM == 2)
            int kk = 0;
            for (int jj = -1; jj <= 1; jj++)  
            for (int ii = -1; ii <= 1; ii++)  
            {
                int index = (jj+1)*3 + (ii+1);
#else
            for (int kk = -1; kk <= 1; kk++)  
            for (int jj = -1; jj <= 1; jj++)  
            for (int ii = -1; ii <= 1; ii++)  
            {
                int index = (kk+1)*9 + (jj+1)*3 + (ii+1);
#endif
                if (nbor(i,j,k,index) == 1)
                {
                    int r = i+ii; int s = j+jj; int t = k+kk;
                    AMREX_D_TERM(cent_hat(i,j,k,0) += (ccent(r,s,t,0) + ii) * vfrac(r,s,t) / nrs(r,s,t);,
                                 cent_hat(i,j,k,1) += (ccent(r,s,t,1) + jj) * vfrac(r,s,t) / nrs(r,s,t);,
                                 cent_hat(i,j,k,2) += (ccent(r,s,t,2) + kk) * vfrac(r,s,t) / nrs(r,s,t););
                }
            }
            AMREX_D_TERM(cent_hat(i,j,k,0) /= nbhd_vol(i,j,k);,
                         cent_hat(i,j,k,1) /= nbhd_vol(i,j,k);,
                         cent_hat(i,j,k,2) /= nbhd_vol(i,j,k););
        } else {
            AMREX_D_TERM(cent_hat(i,j,k,0) = 0.;,
                         cent_hat(i,j,k,1) = 0.;,
                         cent_hat(i,j,k,2) = 0.;);
        }
    });

    // Define Qhat (from Berger and Guliani)
    amrex::ParallelFor(bx, ncomp,  
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (vfrac(i,j,k) > 0.5)
        {
            soln_hat(i,j,k,n) = dUdt_in(i,j,k,n);

        } else if (vfrac(i,j,k) > 0.0) {

#if (AMREX_SPACEDIM == 2)
            int kk = 0;
            for (int jj = -1; jj <= 1; jj++)  
            for (int ii = -1; ii <= 1; ii++)  
            {
                int index = (jj+1)*3 + (ii+1);
#else
            for (int kk = -1; kk <= 1; kk++)  
            for (int jj = -1; jj <= 1; jj++)  
            for (int ii = -1; ii <= 1; ii++)  
            {
                int index = (kk+1)*9 + (jj+1)*3 + (ii+1);
#endif
                if (nbor(i,j,k,index) == 1)
                {
                    int r = i+ii; int s = j+jj; int t = k+kk;
                    soln_hat(i,j,k,n) += dUdt_in(r,s,t,n) * vfrac(r,s,t) / nrs(r,s,t);
                }
            }
            soln_hat(i,j,k,n) /= nbhd_vol(i,j,k);
        } else {
            soln_hat(i,j,k,n) = 1.e40; // NOTE -- we shouldn't end up using this 
        }
    });

    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        dUdt(i,j,k,n) = 0;
    });

    for (int n = 0; n < ncomp; n++)
    {
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (vfrac(i,j,k) > 0.0)
            {
                const auto& slopes_eb = amrex_lim_slopes_eb(i,j,k,n,soln_hat,cent_hat,
                                                            AMREX_D_DECL(fcx,fcy,fcz), flag);
#if (AMREX_SPACEDIM == 2)
            int kk = 0;
            for (int jj = -1; jj <= 1; jj++)  
            for (int ii = -1; ii <= 1; ii++)  
            {
                int index = (jj+1)*3 + (ii+1);
#else
            for (int kk = -1; kk <= 1; kk++)  
            for (int jj = -1; jj <= 1; jj++)  
            for (int ii = -1; ii <= 1; ii++)  
            {
                // Note every cell is at least in its own neighborhood so this will update every (i,j,k)
                int index = (kk+1)*9 + (jj+1)*3 + (ii+1);
#endif
                if (nbor(i,j,k,index) == 1)
                {
                    int r = i+ii; int s = j+jj; int t = k+kk;
                    dUdt(r,s,t,n) += soln_hat(i,j,k,n) + slopes_eb[0] * (ccent(r,s,t,0)-cent_hat(i,j,k,0))
#if (AMREX_SPACEDIM == 3)
                                                       + slopes_eb[2] * (ccent(r,s,t,2)-cent_hat(i,j,k,2))
#endif
                                                       + slopes_eb[1] * (ccent(r,s,t,1)-cent_hat(i,j,k,1));
                } // if nbor(i,j,k,index) == 1
            } // ii,jj,kk loop
            } // if vfrac > 0
        });

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (!flag(i,j,k).isCovered())
            {
                dUdt(i,j,k,n) /= nrs(i,j,k);
            }
        });
    }

    //
    // This tests whether the redistribution procedure was conservative
    //
    { // STRT:SUM OF FINAL DUDT
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
            sum1 += vfrac(i,j,k)*dUdt_in(i,j,k,n);
            sum2 += vfrac(i,j,k)*dUdt(i,j,k,n);
        }
        if (std::abs(sum1-sum2) > 1.e-8 * sum1 && std::abs(sum1-sum2) > 1.e-8)
        {
           printf("SUMS DO NOT MATCH IN STATE REDIST : %f %f ",sum1,sum2);
           amrex::Abort();
        }
      }
    } //  END:SUM OF FINAL DUDT
}
#endif
