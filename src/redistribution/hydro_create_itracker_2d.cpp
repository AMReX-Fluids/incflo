/**
 * \file hydro_create_itracker_2d.cpp
 * \addtogroup Redistribution
 * @{
 *
 */

#include <hydro_redistribution.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 2)

void
Redistribution::MakeITracker ( Box const& bx,
                               Array4<Real const> const& apx,
                               Array4<Real const> const& apy,
                               Array4<Real const> const& vfrac_old,
                               Array4<Real const> const& vfrac_new,
                               Array4<int> const& itracker,
                               Geometry const& lev_geom,
                               Real target_volfrac)
{
#if 0
    int debug_verbose = 0;
#endif

    const Real small_norm_diff = 1.e-8;

    const Box domain = lev_geom.Domain();

    // Note that itracker has 4 components and all are initialized to zero
    // We will add to the first component every time this cell is included in a merged neighborhood,
    //    either by merging or being merged
    // We identify the cells in the remaining three components with the following ordering
    //
    // ^  6 7 8
    // |  4   5
    // j  1 2 3
    //   i --->

    Array<int,9> imap{0,-1,0,1,-1,1,-1,0,1};
    Array<int,9> jmap{0,-1,-1,-1,0,0,1,1,1};

    const auto& is_periodic_x = lev_geom.isPeriodic(0);
    const auto& is_periodic_y = lev_geom.isPeriodic(1);

//  if (debug_verbose > 0)
//      amrex::Print() << " IN MAKE_ITRACKER DOING BOX " << bx << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        itracker(i,j,k,0) = 0;
    });

    Box domain_per_grown = domain;
    if (is_periodic_x) domain_per_grown.grow(0,4);
    if (is_periodic_y) domain_per_grown.grow(1,4);

    Box const& bxg4 = amrex::grow(bx,4);
    Box bx_per_g4= domain_per_grown & bxg4;

    amrex::ParallelFor(bx_per_g4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       if (vfrac_new(i,j,k) > 0.0 && vfrac_new(i,j,k) < target_volfrac)
       {
           Real apnorm, apnorm_inv;
           const Real dapx = apx(i+1,j  ,k  ) - apx(i,j,k);
           const Real dapy = apy(i  ,j+1,k  ) - apy(i,j,k);
           apnorm = std::sqrt(dapx*dapx+dapy*dapy);
           apnorm_inv = 1.0/apnorm;
           Real nx = dapx * apnorm_inv;
           Real ny = dapy * apnorm_inv;

           bool nx_eq_ny = ( (std::abs(nx-ny) < small_norm_diff) ||
                             (std::abs(nx+ny) < small_norm_diff)  ) ? true : false;

           // As a first pass, choose just based on the normal
           if (std::abs(nx) > std::abs(ny))
           {
               if (nx > 0)
                   itracker(i,j,k,1) = 5;
               else
                   itracker(i,j,k,1) = 4;

           } else {
               if (ny > 0)
                   itracker(i,j,k,1) = 7;
               else
                   itracker(i,j,k,1) = 2;
           }

           bool xdir_mns_ok = (is_periodic_x || (i > domain.smallEnd(0)));
           bool xdir_pls_ok = (is_periodic_x || (i < domain.bigEnd(0)  ));
           bool ydir_mns_ok = (is_periodic_y || (j > domain.smallEnd(1)));
           bool ydir_pls_ok = (is_periodic_y || (j < domain.bigEnd(1)  ));

           // Override above logic if trying to reach outside a domain boundary (and non-periodic)
           if ( (!xdir_mns_ok && (itracker(i,j,k,1) == 4)) ||
                (!xdir_pls_ok && (itracker(i,j,k,1) == 5)) )
           {
               itracker(i,j,k,1) = (ny > 0) ? 7 : 2;
           }
           if ( (!ydir_mns_ok && (itracker(i,j,k,1) == 2)) ||
                (!ydir_pls_ok && (itracker(i,j,k,1) == 7)) )
           {
               itracker(i,j,k,1) = (nx > 0) ? 5 : 4;
           }

           // (i,j) merges with at least one cell now
           itracker(i,j,k,0) += 1;

           // (i+ioff,j+joff) is in the nbhd of (i,j)
           int ioff = imap[itracker(i,j,k,1)];
           int joff = jmap[itracker(i,j,k,1)];

           // Sanity check
           if (vfrac_old(i+ioff,j+joff,k) == 0.)
               amrex::Abort(" Trying to merge with covered cell");

           Real sum_vol = vfrac_old(i,j,k) + vfrac_old(i+ioff,j+joff,k);

#if 0
           if (debug_verbose > 0)
               amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vfrac_old(i,j,k) <<
                                 " trying to merge with " << IntVect(i+ioff,j+joff) <<
                                 " with volfrac " << vfrac_old(i+ioff,j+joff,k) <<
                                 " to get new sum_vol " <<  sum_vol << std::endl;
#endif

           // If the merged cell isn't large enough, we try to merge in the other direction
           if (sum_vol < target_volfrac || nx_eq_ny)
           {
               // Original offset was in y-direction, so we will add to the x-direction
               // Note that if we can't because it would go outside the domain, we don't
               if (ioff == 0) {
                   if (nx >= 0 && xdir_pls_ok)
                   {
                       itracker(i,j,k,2) = 5;
                       itracker(i,j,k,0) += 1;
                   }
                   else if (nx <= 0 && xdir_mns_ok)
                   {
                       itracker(i,j,k,2) = 4;
                       itracker(i,j,k,0) += 1;
                   }

               // Original offset was in x-direction, so we will add to the y-direction
               // Note that if we can't because it would go outside the domain, we don't
               } else {
                   if (ny >= 0 && ydir_pls_ok)
                   {
                       itracker(i,j,k,2) = 7;
                       itracker(i,j,k,0) += 1;
                   }
                   else if (ny <= 0 && ydir_mns_ok)
                   {
                       itracker(i,j,k,2) = 2;
                       itracker(i,j,k,0) += 1;
                   }
               }

               if (itracker(i,j,k,0) > 1)
               {
                   // (i+ioff2,j+joff2) is in the nbhd of (i,j)
                   int ioff2 = imap[itracker(i,j,k,2)];
                   int joff2 = jmap[itracker(i,j,k,2)];

                   sum_vol += vfrac_old(i+ioff2,j+joff2,k);
#if 0
                   if (debug_verbose > 0)
                       amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vfrac_old(i,j,k) <<
                                         " trying to ALSO merge with " << IntVect(i+ioff2,j+joff2) <<
                                         " with volfrac " << vfrac_old(i+ioff2,j+joff2,k) <<
                                          " to get new sum_vol " <<  sum_vol << std::endl;
#endif
               }
           }

           // Now we merge in the corner direction if we have already claimed two
           if (itracker(i,j,k,0) == 2)
           {
               // We already have two offsets, and we know they are in different directions
               ioff = imap[itracker(i,j,k,1)] + imap[itracker(i,j,k,2)];
               joff = jmap[itracker(i,j,k,1)] + jmap[itracker(i,j,k,2)];

               if (ioff > 0 && joff > 0)
                   itracker(i,j,k,3) = 8;
               else if (ioff < 0 && joff > 0)
                   itracker(i,j,k,3) = 6;
               else if (ioff > 0 && joff < 0)
                   itracker(i,j,k,3) = 3;
               else
                   itracker(i,j,k,3) = 1;

               // (i,j) merges with at least three cells now
               itracker(i,j,k,0) += 1;

               sum_vol += vfrac_old(i+ioff,j+joff,k);
#if 0
               if (debug_verbose > 0)
                   amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vfrac_old(i,j,k) <<
                                     " trying to ALSO merge with " << IntVect(i+ioff,j+joff) <<
                                     " with volfrac " << vfrac_old(i+ioff,j+joff,k) <<
                                     " to get new sum_vol " <<  sum_vol << std::endl;
#endif
           }
           if (sum_vol < target_volfrac)
           {
#if 0
             amrex::Print() << "Couldnt merge with enough cells to raise volume at " <<
                               IntVect(i,j) << " so stuck with sum_vol " << sum_vol << std::endl;
#endif
             amrex::Abort("Couldnt merge with enough cells to raise volume greater than target_volfrac");
           }
       }
    });
}
#endif
/** @} */
