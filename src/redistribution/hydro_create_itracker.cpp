/*
 * \file hydro_create_itracker.cpp
 * \addtogroup Redistribution
 * @{
 *
 */

#include <hydro_redistribution.H>

#if (AMREX_SPACEDIM == 2)
#include <hydro_create_itracker_2d_K.H>
#else
#include <hydro_create_itracker_3d_K.H>
#endif

using namespace amrex;


void
enforceReciprocity(int i, int j, int k, Array4<int> const& itracker)
{
    auto map = Redistribution::getCellMap();
    // Inverse map
    auto nmap = Redistribution::getInvCellMap();

    // Loop over my neighbors to make sure it's reciprocal, i.e. that my neighbor
    // include me in thier neighborhood too.
    for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
    {
        int ioff = map[0][itracker(i,j,k,i_nbor)];
        int joff = map[1][itracker(i,j,k,i_nbor)];
        int koff = (AMREX_SPACEDIM < 3) ? 0 : map[2][itracker(i,j,k,i_nbor)];

        int ii = i+ioff;
        int jj = j+joff;
        int kk = k+koff;

        if ( Box(itracker).contains(Dim3{ii,jj,kk}) )
        {
            int nbor = itracker(i,j,k,i_nbor);
            int me = nmap[nbor];
            bool found = false;

            // amrex::Print() << "Cell  " << Dim3{i,j,k} << " is (un)covered and merged with neighbor at " << Dim3{i+ioff,j+joff,k+koff} << std::endl;

            // Loop over the neighbor's neighbors to see if I'm already included
            // If not, add me to the neighbor list.
            for (int i_nbor2 = 1; i_nbor2 <= itracker(ii,jj,kk,0); i_nbor2++)
            {
                if ( itracker(ii,jj,kk,i_nbor2) == me ) {
                    // Print()<<IntVect(i,j)<<" is ALREADY A NEIGHBOR!"<<std::endl;
                    found = true;
                    break;
                }
            }
            if ( !found )
            {
                itracker(ii,jj,kk,0) += 1;
                itracker(ii,jj,kk,itracker(ii,jj,kk,0)) = me;
            }
        }
    }
}

void
reverseNbhd(int i, int j, int k, Array4<int> const& itracker)
{
    auto map = Redistribution::getCellMap();
    // Inverse map
    auto nmap = Redistribution::getInvCellMap();

    // Reverse the neighborhood relationship so that now I am in my neighbor's neighborhood
    // and they are no longer in mine
    if (itracker(i,j,k,0) > 1) Abort("reverseNbhd: Nelwy uncovered cell has more than one neighbor");
    for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
    {
        int ioff = map[0][itracker(i,j,k,i_nbor)];
        int joff = map[1][itracker(i,j,k,i_nbor)];
        int koff = (AMREX_SPACEDIM < 3) ? 0 : map[2][itracker(i,j,k,i_nbor)];

        int ii = i+ioff;
        int jj = j+joff;
        int kk = k+koff;

        if ( Box(itracker).contains(Dim3{ii,jj,kk}) )
        {
            int nbor = itracker(i,j,k,i_nbor);
            int me = nmap[nbor];
            bool found = false;

            // amrex::Print() << "Cell  " << Dim3{i,j,k} << " is (un)covered and merged with neighbor at " << Dim3{i+ioff,j+joff,k+koff} << std::endl;

            // Loop over the neighbor's neighbors to see if I'm already included
            // If not, add me to the neighbor list.
            for (int i_nbor2 = 1; i_nbor2 <= itracker(ii,jj,kk,0); i_nbor2++)
            {
                if ( itracker(ii,jj,kk,i_nbor2) == me ) {
                    // Print()<<IntVect(i,j)<<" is ALREADY A NEIGHBOR!"<<std::endl;
                    found = true;
                    Abort("A newly uncovered cell has been included in another neighborhood");
                    break;
                }
            }
            if ( !found )
            {
                itracker(ii,jj,kk,0) += 1;
                itracker(ii,jj,kk,itracker(ii,jj,kk,0)) = me;

                // remove this from my neighborhood.
                itracker(i,j,k,0) -= 1;
                // Newly uncovered cell should now have no neighbors. Otherwise we need to
                // actually remove the neighbor from the array and make adjustments if there's
                // still valid neighbors at higher array indices.
                if (itracker(i,j,k,0) > 0 ) Abort("reverseNbhd: Newly uncovered cell still has neighbors");
            }
        }
    }
}

void
Redistribution::MakeITracker ( Box const& bx,
                               AMREX_D_DECL(Array4<Real const> const& apx_old,
                                            Array4<Real const> const& apy_old,
                                            Array4<Real const> const& apz_old),
                               Array4<Real const> const& vfrac_old,
                               AMREX_D_DECL(Array4<Real const> const& apx_new,
                                            Array4<Real const> const& apy_new,
                                            Array4<Real const> const& apz_new),
                               Array4<Real const> const& vfrac_new,
                               Array4<int> const& itracker,
                               Geometry const& lev_geom,
                               Real target_volfrac,
                               Array4<Real const> const& vel_eb)
{
    int debug_verbose = 0;

    auto map = getCellMap();
    // Inverse map
    auto nmap = getInvCellMap();

    if (debug_verbose > 0)
        amrex::Print() << " IN MAKE_ITRACKER DOING BOX " << bx << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        itracker(i,j,k,0) = 0;
    });

    AMREX_D_TERM(const auto& is_periodic_x = lev_geom.isPeriodic(0);,
                 const auto& is_periodic_y = lev_geom.isPeriodic(1);,
                 const auto& is_periodic_z = lev_geom.isPeriodic(2));
    Box domain_per_grown = lev_geom.Domain();
    AMREX_D_TERM(if (is_periodic_x) domain_per_grown.grow(0,4);,
                 if (is_periodic_y) domain_per_grown.grow(1,4);,
                 if (is_periodic_z) domain_per_grown.grow(2,4));

    Box const& bxg4 = amrex::grow(bx,4);
    Box bx_per_g4= domain_per_grown & bxg4;
    Box domain = lev_geom.Domain();

    amrex::ParallelFor(bx_per_g4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We check for cut-cells in the new geometry
        if ( (vfrac_new(i,j,k) > 0.0 && vfrac_new(i,j,k) < target_volfrac) && vfrac_old(i,j,k) > 0.0)
        {
            normalMerging(i, j, k,
                          AMREX_D_DECL(apx_new, apy_new, apz_new),
                          vfrac_new, itracker,
                          lev_geom, target_volfrac, domain, 
                          AMREX_D_DECL(is_periodic_x, is_periodic_y, is_periodic_z));
        }
        // // WARNING, even with the CFL restriction of MOL, it will NOT always be the case that NU cells
        // // will get a neighbor for default target_volfrac of 0.5. This can happen with a very coarse
        // // simulation due to the single cut cell requirement.
        // else if ( vfrac_new(i,j,k) > 0.0 && vfrac_old(i,j,k) == 0.0)
        // {
        //     // For now, require that newly uncovered cells only have one other cell in it's nbhd
        //     if ( vfrac_new(i,j,k) > target_volfrac )
        //     {
        //         amrex::Warning("WARNING: Newly uncovered cell has volfrac > target_volfrac, suggesting the EB is under-resolved");
        //     }
        //     newlyUncoveredNbhd(i, j, k,
        //                        AMREX_D_DECL(apx_new, apy_new, apz_new),
        //                        vfrac_new, vel_eb, itracker,
        //                        lev_geom, target_volfrac, domain,
        //                        AMREX_D_DECL(is_periodic_x, is_periodic_y, is_periodic_z));
        //     // normalMerging(i, j, k,
        //     //               AMREX_D_DECL(apx_new, apy_new, apz_new),
        //     //               vfrac_new, itracker,
        //     //               lev_geom, target_volfrac, domain,
        //     //               AMREX_D_DECL(is_periodic_x, is_periodic_y, is_periodic_z));
        // }
        else if ( vfrac_old(i,j,k) > 0.0 && vfrac_new(i,j,k) == 0.0)
        {
            // Create a nbhd for cells that become covered...
            // vfrac is only for checking volume of nbhd
            // Probably don't need target_volfrac to match with general case,
            // FIXME? This could result in a NC cell including other NC cells in nbhd depending
            // target_volfrac
            normalMerging(i, j, k,
                          AMREX_D_DECL(apx_old, apy_old, apz_old),
                          vfrac_new, itracker,
                          lev_geom, target_volfrac, domain,
                          AMREX_D_DECL(is_periodic_x, is_periodic_y, is_periodic_z));
        }
    });

// FIXME - need some check to make sure normalMerging doesn't put NU cells into neighborhoods.
    
    // Need a separate loop because this adds to the neighbor's neighborhood
    // probably could alter normalMerging to allow for one loop...
    amrex::ParallelFor(bx_per_g4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // WARNING, even with the CFL restriction of MOL, it will NOT always be the case that NU cells
        // will get a neighbor for default target_volfrac of 0.5. This can happen with a very coarse
        // simulation due to the single cut cell requirement.
        if ( vfrac_new(i,j,k) > 0.0 && vfrac_old(i,j,k) == 0.0)
        {
            // For now, require that newly uncovered cells have one and only one other cell in it's nbhd
            if ( vfrac_new(i,j,k) > target_volfrac )
            {
                amrex::Warning("WARNING: Newly uncovered cell has volfrac > target_volfrac, suggesting the EB is under-resolved");
            }
            newlyUncoveredNbhd(i, j, k,
                               AMREX_D_DECL(apx_new, apy_new, apz_new),
                               vfrac_new, vel_eb, itracker,
                               lev_geom, target_volfrac, domain,
                               AMREX_D_DECL(is_periodic_x, is_periodic_y, is_periodic_z));
            // normalMerging(i, j, k,
            //               AMREX_D_DECL(apx_new, apy_new, apz_new),
            //               vfrac_new, itracker,
            //               lev_geom, target_volfrac, domain,
            //               AMREX_D_DECL(is_periodic_x, is_periodic_y, is_periodic_z));
        }
    });
    
    // // For Newly Uncovered cells, make my neighbor's neighors my own neighbors
    // // since MSRD essentially removes the NU NB cell from the equations with alpha=0
    // amrex::ParallelFor(Box(itracker),
    // [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    // {
    //     if ( (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1. && vfrac_old(i,j,k) == 0.0) ) // Newly uncovered
    //     {
    //         for (int i_nbor = 1; i_nbor <= itr(i,j,k,0); i_nbor++)
    //         {
    //             int ioff = map[0][itracker(i,j,k,i_nbor)];
    //             int joff = map[1][itracker(i,j,k,i_nbor)];
    //             int koff = (AMREX_SPACEDIM < 3) ? 0 : map[2][itracker(i,j,k,i_nbor)];

    //             int ii = i+ioff;
    //             int jj = j+joff;
    //             int kk = k+koff;
    //             if ( Box(itracker).contains(Dim3{ii,jj,kk}) )
    //             {
    //                 if ( itracker(ii, jj, kk, 0) > 0 )
    //                 {
    //                     for (int i_nbor2 = 1; i_nbor2 <= itracker(ii,jj,kk,0); i_nbor2++)
    //                     {
    //                         int ioff2 = map[0][itracker(ii,jj,kk,i_nbor2)];
    //                         int joff2 = map[1][itracker(ii,jj,kk,i_nbor2)];
    //                         int koff2 = (AMREX_SPACEDIM < 3) ? 0 : map[2][itracker(ii,jj,kk,i_nbor2)];

    //                         itracker(i,j,k,0) += 1;
                            
    //                         // WOuls need to extend the mapping to go 2 cells away....
    //                     }
    //             }
    //     }
    // });
    
    // Check uncovered and covered cells, make sure the neighbors also include them.
    // Do this in newlyUncoveredNbhd instead...
//     amrex::ParallelFor(Box(itracker),
//     [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//     {
//         if ( (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1. && vfrac_old(i,j,k) == 0.0 ) // Newly uncovered
//              //|| (vfrac_new(i,j,k) == 0. && vfrac_old(i,j,k) > 0.0) // Newly covered Cells
//          )
//         {
//             int i_nbor = itracker(i,j,k,0);
//             if (i_nbor < 1) Abort("NU cell didn't get a nb");
//             int ii = i+map[0][itracker(i,j,k,i_nbor)];
//             int jj = j+map[1][itracker(i,j,k,i_nbor)];
//             int kk = k+( (AMREX_SPACEDIM < 3) ? 0 : map[2][itracker(i,j,k,i_nbor)] );

//             if ( Box(itracker).contains(Dim3{ii,jj,kk}) )
//             {
//                 // I don;t think this if is needed for enforceReciprocity...
//                 // This was just as a test. To see if when the neighbor took it's own neighbors,
//                 // only then we enforced reciprocity...
//                 // if ( itracker(ii, jj, kk, 0) > 0 )
//                 // {
//                     //enforceReciprocity(i, j, k, itracker);
// // It appears something is off in the math when the neighbor doesn't have it's own neighbors...
//                     reverseNbhd(i, j, k, itracker);
//                 // }
//             }
//         }
//     });

#if 1

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int ii, int jj, int kk) noexcept
    {
        if ( (ii==45 && jj==21 && kk==40) )
        {
            amrex::Print() << "\nInitial Cell Merging" << std::endl;
            
            for ( int i = ii-1; i<= ii+1; i++){
            for ( int j = jj-1; j<= jj+1; j++){
            for ( int k = kk-1; k<= kk+1; k++){ 
            if (itracker(i,j,k) > 0)
            {
            amrex::AllPrint() << "Cell " << Dim3{i,j,k} << " is merged with: ";

            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int ioff = map[0][itracker(i,j,k,i_nbor)];
                int joff = map[1][itracker(i,j,k,i_nbor)];
                int koff = (AMREX_SPACEDIM < 3) ? 0 : map[2][itracker(i,j,k,i_nbor)];

                if (i_nbor > 1)
                {
                    amrex::AllPrint() << ", " << Dim3{i+ioff,j+joff,k+koff};
                } else
                {
                    amrex::AllPrint() << Dim3{i+ioff,j+joff,k+koff};
                }
            }

            amrex::AllPrint() << std::endl;
        }
            }}}
        }
    });
    //amrex::Print() << std::endl;
#endif

#if 0
    amrex::Print() << "Check for all covered cells." << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac_new(i,j,k) == 0. && vfrac_old(i,j,k) > 0.)
        {
            amrex::Print() << "Covered Cell " << Dim3{i,j,k} << std::endl;
        }
    });
    amrex::Print() << std::endl;
#endif

#if 0
    amrex::Print() << "Check for all uncovered cells." << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1. && vfrac_old(i,j,k) == 0.)
        {
            amrex::Print() << "Uncovered Cell " << Dim3{i,j,k} << std::endl;
        }
    });
    amrex::Print() << std::endl;
#endif

#if 0
    amrex::Print() << "Check for all cell that become regular." << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac_old(i,j,k) < 1. && vfrac_new(i,j,k) == 1.)
        {
            amrex::Print() << "New Regular Cell " << Dim3{i,j,k} << std::endl;
        }
    });
    amrex::Print() << std::endl;
#endif

#if 0
    amrex::Print() << "Post Update to Cell Merging" << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        //if(k==8){
        if (itracker(i,j,k) > 0)
        {
            amrex::Print() << "Cell " << Dim3{i,j,k} << " is merged with: ";

            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int ioff = map[0][itracker(i,j,k,i_nbor)];
                int joff = map[1][itracker(i,j,k,i_nbor)];
                int koff = (AMREX_SPACEDIM < 3) ? 0 : map[2][itracker(i,j,k,i_nbor)];

                if (i_nbor > 1)
                {
                    amrex::Print() << ", " << Dim3{i+ioff,j+joff,k+koff};
                } else
                {
                    amrex::Print() << Dim3{i+ioff,j+joff,k+koff};
                }
            }

            amrex::Print() << std::endl;
        }
        //}
    });
    amrex::Print() << std::endl;
#endif

}
/** @} */
