/*
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
			       Array4<EBCellFlag const> const& flag_old,
			       Array4<EBCellFlag const> const& flag_new,
                               Array4<Real const> const& apx_new,
                               Array4<Real const> const& apy_new,
                               Array4<Real const> const& vfrac_old,
			       Array4<Real const> const& vfrac_new,
                               Array4<int> const& itracker,
                               Geometry const& lev_geom,
                               Real target_volfrac)
{
#if 1
    int debug_verbose = 1;
#endif

//    constexpr Real eps = 1.e-8; 
    constexpr Real small_norm_diff = 1.e-8;

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
    // Note the first component should never get used (it refers to the cell itself, the blank in the center above)
    Array<int,9> imap{0,-1,0,1,-1,1,-1,0,1};
    Array<int,9> jmap{0,-1,-1,-1,0,0,1,1,1};
    // CEG for reciprocity with newly uncovered cells
    //Array<int,9> nmap{0,0,7,0,5,4,0,2,0};

    const auto& is_periodic_x = lev_geom.isPeriodic(0);
    const auto& is_periodic_y = lev_geom.isPeriodic(1);

    if (debug_verbose > 0)
        amrex::Print() << " IN MAKE_ITRACKER DOING BOX " << bx << std::endl;

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
	// We check for cut-cells in the new geometry 
	if (!flag_new(i,j,k).isCovered()  && (vfrac_new(i,j,k) < target_volfrac || vfrac_old(i,j,k) < target_volfrac) )
	{
// 	   //
// 	   // Choose based on EB normal.
// 	   //
//            Real apnorm, apnorm_inv, dapx, dapy;
           
//            dapx = apx_new(i+1,j,k) - apx_new(i,j,k);
//            dapy = apy_new(i,j+1,k) - apy_new(i,j,k);

// 	   //Real deltaV = vfrac_new(i,j,k) - vfrac_old(i,j,k);
	   
//            // // CEG: Matthew's first pass at trying to select for cut cells behind the flow
//            // // of the EB, particularly relevant to covered-to-cut cells that require reciprocity
// 	   // // Start with normal and flip sign to select movement of EB
//            // if ( (vfrac_new(i,j,k) - vfrac_old(i,j,k) > 0.0) )
//            // {
//            //      dapx = apx_new(i,j,k) - apx_new(i+1, j, k);
//            //      dapy = apy_new(i,j,k) - apy_new(i,j+1, k);
//            // } else {
//            //      dapx = apx_new(i+1,j,k) - apx_new(i,j,k);
//            //      dapy = apy_new(i,j+1,k) - apy_new(i,j,k);
//            // }
 
//            apnorm = std::sqrt(dapx*dapx+dapy*dapy);
//            apnorm_inv = 1.0/apnorm;
//            Real nx = dapx * apnorm_inv;
//            Real ny = dapy * apnorm_inv;

//            bool nx_eq_ny = ( (std::abs(nx-ny) < small_norm_diff) ||
//                              (std::abs(nx+ny) < small_norm_diff)  ) ? true : false;

//            // Original: As a first pass, choose just based on the normal 
//            if (std::abs(nx) > std::abs(ny))
//            {
//                if (nx > 0)
//                    itracker(i,j,k,1) = 5;
//                else
//                    itracker(i,j,k,1) = 4;

//            } else {
//                if (ny > 0)
//                    itracker(i,j,k,1) = 7;
//                else
//                    itracker(i,j,k,1) = 2;
//            }
	   
// 	   // Select with bias for EB movement
// 	   // // need a tolerance/eps here...
//            // if (std::abs(nx) > std::abs(ny)) {
//            //     if (nx > 0) {
//            //         if (vfrac_new(i-1,j,k) < eps && vfrac_old(i-1,j,k) < eps){
//            //             itracker(i,j,k,1) = 5;
//            //         } else {
//            //             itracker(i,j,k,1) = 4;
//            //         }
//            //     } else {
//            //         if (vfrac_new(i+1,j,k) < eps && vfrac_old(i+1,j,k) < eps){
//            //             itracker(i,j,k,1) = 4;
//            //         } else {
//            //             itracker(i,j,k,1) = 5;
//            //         }
//            //     }
//            // } else {
//            //     if (ny > 0) {
//            //         if (vfrac_new(i,j-1,k) < eps && vfrac_old(i,j-1,k) < eps){
//            //             itracker(i,j,k,1) = 7;
//            //         } else {
//            //             itracker(i,j,k,1) = 2;
//            //         }
//            //     } else {
//            //         if (vfrac_new(i,j+1,k) < eps && vfrac_old(i,j+1,k) < eps){
//            //             itracker(i,j,k,1) = 2; 
//            //         } else {
//            //             itracker(i,j,k,1) = 7;
//            //         }
//            //     }
//            // }

//            bool xdir_mns_ok = (is_periodic_x || (i > domain.smallEnd(0)));
//            bool xdir_pls_ok = (is_periodic_x || (i < domain.bigEnd(0)  ));
//            bool ydir_mns_ok = (is_periodic_y || (j > domain.smallEnd(1)));
//            bool ydir_pls_ok = (is_periodic_y || (j < domain.bigEnd(1)  ));

//            // Override above logic if trying to reach outside a domain boundary (and non-periodic)
//            if ( (!xdir_mns_ok && (itracker(i,j,k,1) == 4)) ||
//                 (!xdir_pls_ok && (itracker(i,j,k,1) == 5)) )
//            {
//                itracker(i,j,k,1) = (ny > 0) ? 7 : 2;
//            }
//            if ( (!ydir_mns_ok && (itracker(i,j,k,1) == 2)) ||
//                 (!ydir_pls_ok && (itracker(i,j,k,1) == 7)) )
//            {
//                itracker(i,j,k,1) = (nx > 0) ? 5 : 4;
//            }

//            // (i,j) merges with at least one cell now
//            itracker(i,j,k,0) += 1;

//            // (i+ioff,j+joff) is in the nbhd of (i,j)
//            int ioff = imap[itracker(i,j,k,1)];
//            int joff = jmap[itracker(i,j,k,1)];

// #if 1
//            Real sum_vol = vfrac_new(i,j,k) + vfrac_new(i+ioff,j+joff,k);
           
//            if ( i==16 && j==8 )//( debug_verbose > 0 )
//               amrex::Print() << "Cell " << IntVect(i,j) << " with vfrac " << vfrac_new(i,j,k) <<
//                                 " merge " << IntVect(i+ioff,j+joff) <<
//                                 " with vfrac " << vfrac_new(i+ioff,j+joff,k) <<
//                                 " to get new vfrac " <<  sum_vol << std::endl;
// #endif

                 
// #if 1
//            // If the merged cell isn't large enough, we try to merge in the other direction
//            if (sum_vol < target_volfrac || nx_eq_ny)
//            {
//                // Original offset was in y-direction, so we will add to the x-direction
//                // Note that if we can't because it would go outside the domain, we don't
//                if (ioff == 0) {
//                    if (nx >= 0 && xdir_pls_ok)
//                    {
//                        itracker(i,j,k,2) = 5;
//                        itracker(i,j,k,0) += 1;
//                    }
//                    else if (nx <= 0 && xdir_mns_ok)
//                    {
//                        itracker(i,j,k,2) = 4;
//                        itracker(i,j,k,0) += 1;
//                    }

//                // Original offset was in x-direction, so we will add to the y-direction
//                // Note that if we can't because it would go outside the domain, we don't
//                } else {
//                    if (ny >= 0 && ydir_pls_ok)
//                    {
//                        itracker(i,j,k,2) = 7;
//                        itracker(i,j,k,0) += 1;
//                    }
//                    else if (ny <= 0 && ydir_mns_ok)
//                    {
//                        itracker(i,j,k,2) = 2;
//                        itracker(i,j,k,0) += 1;
//                    }
//                }

//                if (itracker(i,j,k,0) > 1)
//                {
//                    // (i+ioff2,j+joff2) is in the nbhd of (i,j)
//                    int ioff2 = imap[itracker(i,j,k,2)];
//                    int joff2 = jmap[itracker(i,j,k,2)];

//                    sum_vol += vfrac_new(i+ioff2,j+joff2,k);
// #if 1
//                    if ( i==16 && j==8 ) // (debug_verbose > 0 )
//                        amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vfrac_new(i,j,k) <<
//                                          " trying to ALSO merge with " << IntVect(i+ioff2,j+joff2) <<
//                                          " with volfrac " << vfrac_new(i+ioff2,j+joff2,k) <<
//                                           " to get new sum_vol " <<  sum_vol << std::endl;
// #endif
//                }
//            }
// #endif
 
//            // Now we merge in the corner direction if we have already claimed two
//            if (itracker(i,j,k,0) == 2)
//            {
//                // We already have two offsets, and we know they are in different directions
// 	       // don't shadow
//                ioff = imap[itracker(i,j,k,1)] + imap[itracker(i,j,k,2)];
//                joff = jmap[itracker(i,j,k,1)] + jmap[itracker(i,j,k,2)];
               
//                if (ioff > 0 && joff > 0)
//                    itracker(i,j,k,3) = 8;
//                else if (ioff < 0 && joff > 0)
//                    itracker(i,j,k,3) = 6;
//                else if (ioff > 0 && joff < 0)
//                    itracker(i,j,k,3) = 3;
//                else
//                    itracker(i,j,k,3) = 1;

//                // (i,j) merges with at least three cells now
//                itracker(i,j,k,0) += 1;

//                sum_vol += vfrac_new(i+ioff,j+joff,k);
// #if 1
//                if ( i==16 && j==8 ) // (debug_verbose > 0 )
//                    amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vfrac_new(i,j,k) <<
//                                      " trying to ALSO merge with " << IntVect(i+ioff,j+joff) <<
//                                      " with volfrac " << vfrac_new(i+ioff,j+joff,k) <<
//                                      " to get new sum_vol " <<  sum_vol << std::endl;
// #endif
//            }

           //
	   // Central Merging
	   //
	    Real sum_vol = vfrac_new(i,j,k);
	    
	    //int label[3][3] = {{1,2,3}, {4,0,5}, {6,7,8}};
	    int label[3][3] = {{1,4,6}, {2,0,7}, {3,5,8}};	   
	   int kk=0;
	   for(int jj(-1); jj<=1; jj++) {
	       for(int ii(-1); ii<=1; ii++) {
		   if (flag_new(i,j,k).isConnected(ii,jj,kk))
		   {
		       if ( !(ii==0 && jj==0)
			   &&(is_periodic_x || (domain.smallEnd(0)<=i+ii && i+ii<=domain.bigEnd(0)))
			   &&(is_periodic_y || (domain.smallEnd(1)<=j+jj && j+jj<=domain.bigEnd(1))))
		       {
			   // take all valid interior cells in 3x3 neighborhood centered
			   // excluding i,j itself
			   // this includes any newly uncovered cells
			   itracker(i,j,k,0) += 1;
			   itracker(i,j,k,itracker(i,j,k,0)) = label[ii+1][jj+1];
			   sum_vol += vfrac_new(i+ii,j+jj,k);

			   //FIXME -- still need to enfore reciprocity for newly uncovered cells
//fixme
			   if ( i==16 && j==8 ){
			       Print()<<"Including cell ("<<ii<<","<<jj<<"). label: "
				      <<label[ii+1][jj+1]<<std::endl;
			   }
		       }
		   }
		   //FIXME - do we really want to try to include newly covered cells here?
		   // they can't fully participate in SRD because we don't give them a nbhd...
		   else if (flag_old(i,j,k).isConnected(ii,jj,kk))
		   {
                       // Add newly covered cells to the neighborhood
		       if (  (is_periodic_x || (domain.smallEnd(0)<=i+ii && i+ii<=domain.bigEnd(0)))
			   &&(is_periodic_y || (domain.smallEnd(1)<=j+jj && j+jj<=domain.bigEnd(1))))
		       {
			   itracker(i,j,k,0) += 1;
			   itracker(i,j,k,itracker(i,j,k,0)) = label[ii+1][jj+1];
//fixme
			   if ( i==16 && j==8 ){
			       Print()<<"Including newly covered cell ("<<ii<<","<<jj<<"). label: "
				      <<label[ii+1][jj+1]<<std::endl;
			   }
		       }		       
		   }
	       }
	   }
#if 1
           if (sum_vol < target_volfrac)
           {
             amrex::Print() << "Couldnt merge with enough cells to raise volume at " <<
                               IntVect(i,j) << " so stuck with sum_vol " << sum_vol << std::endl;
             amrex::Warning("Couldnt merge with enough cells to raise volume greater than target_volfrac");
	     if (sum_vol < Real(0.5))
	     {
		 amrex::Abort("Couldnt merge with enough cells to raise volume greater than 0.5");
	     }
           }
#endif
       }
    });


#if 1
    amrex::Print() << "\nInitial Cell Merging" << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
//fixme
	if ( i==16 && j==8 )
        //if (itracker(i,j,k) > 0)
        {
            amrex::Print() << "Cell " << IntVect(i,j) << " is merged with "
			   <<itracker(i,j,k,0)<<" cells :\n"; 

            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int ioff = imap[itracker(i,j,k,i_nbor)];
                int joff = jmap[itracker(i,j,k,i_nbor)];  
                
                if (i_nbor > 1)
                {
                    amrex::Print() << ", " << IntVect(i+ioff, j+joff);
                } else
                {
                    amrex::Print() << IntVect(i+ioff, j+joff);
                }
            }

            amrex::Print() << std::endl;
        }
    });    
    amrex::Print() << std::endl;
#endif

    // //
    // // Is this meant to implement neighborhood reciprocity??? 
    // // Central merging takes care of some reciprocity by virture of included all uncovered cells
    // //
    // // FIXME - need to make sure we don't reach outside domain here...
    // // I think this should happen before including the middle cell...
    // amrex::ParallelFor(Box(itracker),
    // [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    // {
    // 	// Check if my neighborhood includes any covered cells. If yes, add this
    //  // cell to the covered cells nbhd. 
    //     if (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < target_volfrac)
    //     {
    //         for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
    //         {
    // FIXME, doesn't this need i_nbor instead of hard-coded 1???
    //             int ioff = imap[itracker(i,j,k,1)];
    //             int joff = jmap[itracker(i,j,k,1)];   
           
    //             if (vfrac_new(i+ioff,j+joff,k) < eps )
    //             {
    //                 amrex::Print() << "Cell  " << IntVect(i,j) << " is merged with a covered neighbor at " << IntVect(i+ioff,j+joff) << std::endl;
    //                 itracker(i+ioff,j+joff,k,0) += 1;
    //                 itracker(i+ioff,j+joff,k,itracker(i+ioff,j+joff,k,0)) = nmap[itracker(i,j,k,1)];
    //             }
    //         }
    //     }
    // 	// Add neighboors to uncovered cells
    //     if (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1. && vfrac_old(i,j,k) < eps)
    //     {
    //         for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)            
    //         {
        // FIXME, doesn't this need i_nbor instead of hard-coded 1???
    // 		amrex::Print() << "Uncovered Cell " << IntVect(i,j) << std::endl;   
    //             int ioff = imap[itracker(i,j,k,1)];
    //             int joff = jmap[itracker(i,j,k,1)];   
           
    // 		itracker(i+ioff,j+joff,k,0) += 1;
    // 		itracker(i+ioff,j+joff,k,itracker(i+ioff,j+joff,k,0)) = nmap[itracker(i,j,k,1)];
    //         }
    //     }

    // });

#if 1
    amrex::Print() << "Check for all newly covered cells." << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
	if (flag_new(i,j,k).isCovered() && !flag_old(i,j,k).isCovered())
        {
            amrex::Print() << "Covered Cell " << IntVect(i,j) << std::endl;   
        }
    });    
    amrex::Print() << std::endl;
#endif

#if 1
    amrex::Print() << "Check for all newly uncovered cells." << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
	if (!flag_new(i,j,k).isCovered() && flag_old(i,j,k).isCovered())
        {
            amrex::Print() << "Uncovered Cell " << IntVect(i,j) << std::endl;   
        }
    });    
    amrex::Print() << std::endl;
#endif

#if 1
    amrex::Print() << "Check for all cell that become regular." << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (!flag_old(i,j,k).isRegular() && flag_new(i,j,k).isRegular())
        {
            amrex::Print() << "New Regular Cell " << IntVect(i,j) << std::endl;   
        }
    });    
    amrex::Print() << std::endl;
#endif

#if 0
    amrex::Print() << "Post Update to Cell Merging" << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (itracker(i,j,k) > 0)
        {
            amrex::Print() << "Cell " << IntVect(i,j) << " is merged with: ";   

            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int ioff = imap[itracker(i,j,k,i_nbor)];
                int joff = jmap[itracker(i,j,k,i_nbor)];  
                
                if (i_nbor > 1)
                {
                    amrex::Print() << ", " << IntVect(i+ioff, j+joff);
                } else
                {
                    amrex::Print() << IntVect(i+ioff, j+joff);
                }
            }

            amrex::Print() << std::endl;
        }
    });    
    amrex::Print() << std::endl;
#endif

}
#endif
/** @} */
