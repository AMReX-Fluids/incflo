/**
 * \file hydro_create_itracker_3d.cpp
 * \addtogroup Redistribution
 * @{
 *
 */

#include <hydro_redistribution.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 3)

amrex::Array<amrex::Array<int,27>,AMREX_SPACEDIM>
Redistribution::getCellMap()
{
    //
    // Note that itracker has XX components and all are initialized to zero
    // We will add to the first component every time this cell is included in a merged neighborhood,
    //    either by merging or being merged
    // We identify the cells in the remaining three components with the following ordering
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
    Array<int,27>    imap{0,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    Array<int,27>    jmap{0,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    Array<int,27>    kmap{0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    amrex::Array<amrex::Array<int,27>,AMREX_SPACEDIM> map{imap,jmap,kmap};
    return map;
}

amrex::Array<int,27>
Redistribution::getInvCellMap()
{
    //
    // Create the inverse cell map for
    //
    //    at k-1   |   at k  |   at k+1
    //
    // ^  15 16 17 |  6 7 8  |  24 25 26
    // |  12 13 14 |  4   5  |  21 22 23
    // j  9  10 11 |  1 2 3  |  18 19 20
    //   i --->
    //
    Array<int,27> nmap{0,8,7,6,5,4,3,2,1,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9};
    return nmap;
}

void
Redistribution::normalMerging ( int i, int j, int k,
                                Array4<Real const> const& apx,
                                Array4<Real const> const& apy,
                                Array4<Real const> const& apz,
                                Array4<Real const> const& vfrac,
                                Array4<int> const& itracker,
                                Geometry const& geom,
                                Real target_volfrac)
{
#if 0
    bool debug_print = false;
#endif

    using defaults::small_norm_diff;

    const Box domain = geom.Domain();
    const auto& is_periodic_x = geom.isPeriodic(0);
    const auto& is_periodic_y = geom.isPeriodic(1);
    const auto& is_periodic_z = geom.isPeriodic(2);

    Real apnorm, apnorm_inv;
    const Real dapx = apx(i+1,j  ,k  ) - apx(i,j,k);
    const Real dapy = apy(i  ,j+1,k  ) - apy(i,j,k);
    const Real dapz = apz(i  ,j  ,k+1) - apz(i,j,k);
    apnorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
    apnorm_inv = 1.0/apnorm;
    Real nx = dapx * apnorm_inv;
    Real ny = dapy * apnorm_inv;
    Real nz = dapz * apnorm_inv;

    bool nx_eq_ny = ( (std::abs(nx-ny) < small_norm_diff) ||
                      (std::abs(nx+ny) < small_norm_diff)  ) ? true : false;
    bool nx_eq_nz = ( (std::abs(nx-nz) < small_norm_diff) ||
                      (std::abs(nx+nz) < small_norm_diff)  ) ? true : false;
    bool ny_eq_nz = ( (std::abs(ny-nz) < small_norm_diff) ||
                      (std::abs(ny+nz) < small_norm_diff)  ) ? true : false;

    bool xdir_mns_ok = (is_periodic_x || (i > domain.smallEnd(0)));
    bool xdir_pls_ok = (is_periodic_x || (i < domain.bigEnd(0)  ));
    bool ydir_mns_ok = (is_periodic_y || (j > domain.smallEnd(1)));
    bool ydir_pls_ok = (is_periodic_y || (j < domain.bigEnd(1)  ));
    bool zdir_mns_ok = (is_periodic_z || (k > domain.smallEnd(2)));
    bool zdir_pls_ok = (is_periodic_z || (k < domain.bigEnd(2)  ));

    // x-component of normal is greatest
    if ( (std::abs(nx) > std::abs(ny)) &&
         (std::abs(nx) > std::abs(nz)) )
    {
        if (nx > 0)
            itracker(i,j,k,1) = 5;
        else
            itracker(i,j,k,1) = 4;

        // y-component of normal is greatest
    } else if ( (std::abs(ny) >= std::abs(nx)) &&
                (std::abs(ny) > std::abs(nz)) )
    {
        if (ny > 0)
            itracker(i,j,k,1) = 7;
        else

            itracker(i,j,k,1) = 2;
        // z-component of normal is greatest
    } else {
        if (nz > 0)
            itracker(i,j,k,1) = 22;
        else
            itracker(i,j,k,1) = 13;
    }

    // Override above logic if trying to reach outside a domain boundary (and non-periodic)
    if ( (!xdir_mns_ok && (itracker(i,j,k,1) == 4)) ||
         (!xdir_pls_ok && (itracker(i,j,k,1) == 5)) )
    {
        if ( (std::abs(ny) > std::abs(nz)) )
            itracker(i,j,k,1) = (ny > 0) ? 7 : 2;
        else
            itracker(i,j,k,1) = (nz > 0) ? 22 : 13;
    }

    if ( (!ydir_mns_ok && (itracker(i,j,k,1) == 2)) ||
         (!ydir_pls_ok && (itracker(i,j,k,1) == 7)) )
    {
        if ( (std::abs(nx) > std::abs(nz)) )
            itracker(i,j,k,1) = (nx > 0) ? 5 : 4;
        else
            itracker(i,j,k,1) = (nz > 0) ? 22 : 13;
    }

    if ( (!zdir_mns_ok && (itracker(i,j,k,1) == 13)) ||
         (!zdir_pls_ok && (itracker(i,j,k,1) == 22)) )
    {
        if ( (std::abs(nx) > std::abs(ny)) )
            itracker(i,j,k,1) = (nx > 0) ? 5 : 4;
        else
            itracker(i,j,k,1) = (ny > 0) ? 7 : 2;
    }

    // (i,j,k) merges with at least one cell now
    itracker(i,j,k,0) += 1;

    // (i+ioff,j+joff,k+koff) is now the first cell in the nbhd of (i,j,k)
    auto map = getCellMap();
    int ioff = map[0][itracker(i,j,k,1)];
    int joff = map[1][itracker(i,j,k,1)];
    int koff = map[2][itracker(i,j,k,1)];

    // Sanity check
    if (vfrac(i+ioff,j+joff,k+koff) == 0.)
    {
        amrex::Print() << "Cell " << IntVect(i,j,k) << " is trying to merge with cell " << IntVect(i+ioff,j+joff,k+koff) << std::endl;
        amrex::Abort("normalMerging(): Trying to merge with covered cell");
    }

    Real sum_vol = vfrac(i,j,k) + vfrac(i+ioff,j+joff,k+koff);

#if 0
    if (debug_print)
        amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vfrac(i,j,k) <<
            " trying to merge with " << IntVect(i+ioff,j+joff,k+koff) <<
            " with volfrac " << vfrac(i+ioff,j+joff,k+koff) <<
            " to get new sum_vol " <<  sum_vol << std::endl;
#endif

    bool just_broke_symmetry = ( ( (joff == 0 && koff == 0) && (nx_eq_ny || nx_eq_nz) ) ||
                                 ( (ioff == 0 && koff == 0) && (nx_eq_ny || ny_eq_nz) ) ||
                                 ( (ioff == 0 && joff == 0) && (nx_eq_nz || ny_eq_nz) ) );

    // If the merged cell isn't large enough, or if we broke symmetry by the current merge,
    // we merge in one of the other directions.  Note that the direction of the next merge
    // is first set by a symmetry break, but if that isn't happening, we choose the next largest normal

    if ( (sum_vol < target_volfrac) || just_broke_symmetry )
    {
        // Original offset was in x-direction
        if (joff == 0 && koff == 0)
        {
            if (nx_eq_ny) {
                itracker(i,j,k,2) = (ny > 0) ? 7 : 2;
            } else if (nx_eq_nz) {
                itracker(i,j,k,2) = (nz > 0) ? 22 : 13;
            } else if ( (std::abs(ny) > std::abs(nz)) ) {
                itracker(i,j,k,2) = (ny > 0) ? 7 : 2;
            } else {
                itracker(i,j,k,2) = (nz > 0) ? 22 : 13;
            }

            // Original offset was in y-direction
        } else if (ioff == 0 && koff == 0)
        {
            if (nx_eq_ny) {
                itracker(i,j,k,2) = (nx > 0) ? 5 : 4;
            } else if (ny_eq_nz) {
                itracker(i,j,k,2) = (nz > 0) ? 22 : 13;
            } else if ( (std::abs(nx) > std::abs(nz)) ) {
                itracker(i,j,k,2) = (nx > 0) ? 5 : 4;
            } else {
                itracker(i,j,k,2) = (nz > 0) ? 22 : 13;
            }

            // Original offset was in z-direction
        } else if (ioff == 0 && joff == 0)
        {
            if (nx_eq_nz) {
                itracker(i,j,k,2) = (nx > 0) ? 5 : 4;
            } else if (ny_eq_nz) {
                itracker(i,j,k,2) = (ny > 0) ? 7 : 2;
            } else if ( (std::abs(nx) > std::abs(ny)) ) {
                itracker(i,j,k,2) = (nx > 0) ? 5 : 4;
            } else {
                itracker(i,j,k,2) = (ny > 0) ? 7 : 2;
            }
        }

        // (i,j,k) merges with at least two cells now
        itracker(i,j,k,0) += 1;

        // (i+ioff2,j+joff2,k+koff2) is in the nbhd of (i,j,k)
        int ioff2 = map[0][itracker(i,j,k,2)];
        int joff2 = map[1][itracker(i,j,k,2)];
        int koff2 = map[2][itracker(i,j,k,2)];

        sum_vol += vfrac(i+ioff2,j+joff2,k+koff2);
#if 0
        if (debug_print)
            amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vfrac(i,j,k) <<
                " trying to ALSO merge with " << IntVect(i+ioff2,j+joff2,k+koff2) <<
                " with volfrac " << vfrac(i+ioff2,j+joff2,k+koff2) <<
                " to get new sum_vol " <<  sum_vol << std::endl;
#endif
    }

    // If the merged cell has merged in two directions, we now merge in the corner direction within the current plane
    if (itracker(i,j,k,0) >= 2)
    {
        // We already have two offsets, and we know they are in different directions
        ioff = map[0][itracker(i,j,k,1)] + map[0][itracker(i,j,k,2)];
        joff = map[1][itracker(i,j,k,1)] + map[1][itracker(i,j,k,2)];
        koff = map[2][itracker(i,j,k,1)] + map[2][itracker(i,j,k,2)];

        // Both nbors are in the koff=0 plane
        if (koff == 0)
        {
            if (ioff > 0 && joff > 0)
                itracker(i,j,k,3) = 8;
            else if (ioff < 0 && joff > 0)
                itracker(i,j,k,3) = 6;
            else if (ioff > 0 && joff < 0)
                itracker(i,j,k,3) = 3;
            else
                itracker(i,j,k,3) = 1;

            // Both nbors are in the joff=0 plane
        } else if (joff == 0) {
            if (ioff > 0 && koff > 0)
                itracker(i,j,k,3) = 23;
            else if (ioff < 0 && koff > 0)
                itracker(i,j,k,3) = 21;
            else if (ioff > 0 && koff < 0)
                itracker(i,j,k,3) = 14;
            else
                itracker(i,j,k,3) = 12;

            // Both nbors are in the ioff=0 plane
        } else {
            if (joff > 0 && koff > 0)
                itracker(i,j,k,3) = 25;
            else if (joff < 0 && koff > 0)
                itracker(i,j,k,3) = 19;
            else if (joff > 0 && koff < 0)
                itracker(i,j,k,3) = 16;
            else
                itracker(i,j,k,3) = 10;
        }

        // (i,j,k) merges with at least three cells now
        itracker(i,j,k,0) += 1;

        sum_vol += vfrac(i+ioff,j+joff,k+koff);

#if 0
        int ioff3 = map[0][itracker(i,j,k,3)];
        int joff3 = map[1][itracker(i,j,k,3)];
        int koff3 = map[2][itracker(i,j,k,3)];
        if (debug_print)
            amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vfrac(i,j,k) <<
                " trying to ALSO merge with " << IntVect(i+ioff3,j+joff3,k+koff3) <<
                " with volfrac " << vfrac(i+ioff3,j+joff3,k+koff3) << std::endl;
#endif

        // All nbors are currently in one of three planes
        just_broke_symmetry = ( ( (koff == 0) && (nx_eq_nz || ny_eq_nz) ) ||
                                ( (joff == 0) && (nx_eq_ny || ny_eq_nz) ) ||
                                ( (ioff == 0) && (nx_eq_ny || nx_eq_nz) ) );

        // If with a nbhd of four cells we have still not reached vfrac > target_volfrac, we add another four
        //    cells to the nbhd to make a 2x2x2 block.  We use the direction of the remaining
        //    normal to know whether to go lo or hi in the new direction.
        if (sum_vol < target_volfrac || just_broke_symmetry)
        {
#if 0
            if (debug_print)
                if (just_broke_symmetry)
                    amrex::Print() << "Expanding neighborhood of " << IntVect(i,j,k) <<
                        " from 4 to 8 since we just broke symmetry with the last merge " << std::endl;
                else
                    amrex::Print() << "Expanding neighborhood of " << IntVect(i,j,k) <<
                        " from 4 to 8 since sum_vol with 4 was only " << sum_vol << " " << std::endl;
#endif
            // All nbors are currently in the koff=0 plane
            if (koff == 0)
            {
                if (nz > 0)
                {
                    itracker(i,j,k,4) = 22;

                    if (ioff > 0)
                        itracker(i,j,k,5) =  23;
                    else
                        itracker(i,j,k,5) =  21;
                    if (joff > 0)
                        itracker(i,j,k,6) =  25;
                    else
                        itracker(i,j,k,6) =  19;

                    if (ioff > 0 && joff > 0) {
                        itracker(i,j,k,7) = 26;
                    } else if (ioff < 0 && joff > 0) {
                        itracker(i,j,k,7) = 24;
                    } else if (ioff > 0 && joff < 0) {
                        itracker(i,j,k,7) = 20;
                    } else {
                        itracker(i,j,k,7) = 18;
                    }
                } else { // nz <= 0

                    itracker(i,j,k,4) = 13;

                    if (ioff > 0)
                        itracker(i,j,k,5) =  14;
                    else
                        itracker(i,j,k,5) =  12;
                    if (joff > 0)
                        itracker(i,j,k,6) =  16;
                    else
                        itracker(i,j,k,6) =  10;

                    if (ioff > 0 && joff > 0) {
                        itracker(i,j,k,7) = 17;
                    } else if (ioff < 0 && joff > 0) {
                        itracker(i,j,k,7) = 15;
                    } else if (ioff > 0 && joff < 0) {
                        itracker(i,j,k,7) = 11;
                    } else {
                        itracker(i,j,k,7) =  9;
                    }
                }
            } else if (joff == 0) {
                if (ny > 0)
                {
                    itracker(i,j,k,4) = 7;

                    if (ioff > 0)
                        itracker(i,j,k,5) =  8;
                    else
                        itracker(i,j,k,5) =  6;
                    if (koff > 0)
                        itracker(i,j,k,6) =  25;
                    else
                        itracker(i,j,k,6) =  16;

                    if (ioff > 0 && koff > 0) {
                        itracker(i,j,k,7) = 26;
                    } else if (ioff < 0 && koff > 0) {
                        itracker(i,j,k,7) = 24;
                    } else if (ioff > 0 && koff < 0) {
                        itracker(i,j,k,7) = 17;
                    } else {
                        itracker(i,j,k,7) = 15;
                    }

                } else { // ny <= 0

                    itracker(i,j,k,4) = 2;

                    if (ioff > 0)
                        itracker(i,j,k,5) =  3;
                    else
                        itracker(i,j,k,5) =  1;
                    if (koff > 0)
                        itracker(i,j,k,6) =  19;
                    else
                        itracker(i,j,k,6) =  10;

                    if (ioff > 0 && koff > 0) {
                        itracker(i,j,k,7) = 20;
                    } else if (ioff < 0 && koff > 0) {
                        itracker(i,j,k,7) = 18;
                    } else if (ioff > 0 && koff < 0) {
                        itracker(i,j,k,7) = 11;
                    } else {
                        itracker(i,j,k,7) =  9;
                    }
                }
            } else if (ioff == 0) {

                if (nx > 0)
                {
                    itracker(i,j,k,4) = 5;

                    if (joff > 0)
                        itracker(i,j,k,5) =  8;
                    else
                        itracker(i,j,k,5) =  3;
                    if (koff > 0)
                        itracker(i,j,k,6) =  23;
                    else
                        itracker(i,j,k,6) =  14;

                    if (joff > 0 && koff > 0) {
                        itracker(i,j,k,7) = 26;
                    } else if (joff < 0 && koff > 0) {
                        itracker(i,j,k,7) = 20;
                    } else if (joff > 0 && koff < 0) {
                        itracker(i,j,k,7) = 17;
                    } else {
                        itracker(i,j,k,7) = 11;
                    }
                } else { // nx <= 0

                    itracker(i,j,k,4) = 4;

                    if (joff > 0)
                        itracker(i,j,k,5) =  6;
                    else
                        itracker(i,j,k,5) =  1;
                    if (koff > 0)
                        itracker(i,j,k,6) =  21;
                    else
                        itracker(i,j,k,6) =  12;

                    if (joff > 0 && koff > 0) {
                        itracker(i,j,k,7) = 24;
                    } else if (joff < 0 && koff > 0) {
                        itracker(i,j,k,7) = 18;
                    } else if (joff > 0 && koff < 0) {
                        itracker(i,j,k,7) = 15;
                    } else {
                        itracker(i,j,k,7) =  9;
                    }
                }
            }

            for(int n(4); n<8; n++){
                int ioffn = map[0][itracker(i,j,k,n)];
                int joffn = map[1][itracker(i,j,k,n)];
                int koffn = map[2][itracker(i,j,k,n)];
                sum_vol += vfrac(i+ioffn,j+joffn,k+koffn);
#if 0
                if (debug_print)
                    amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vfrac(i,j,k) <<
                        " trying to ALSO merge with " << IntVect(i+ioffn,j+joffn,k+koffn) <<
                        " with volfrac " << vfrac(i+ioffn,j+joffn,k+koffn) <<
                        " to get new sum_vol " <<  sum_vol << std::endl;
#endif
            }

            // (i,j,k) has a 2x2x2 neighborhood now
            itracker(i,j,k,0) += 4;
        }
    }
    if (sum_vol < target_volfrac)
    {
#if 0
        amrex::Print() << "Couldnt merge with enough cells to raise volume at " <<
            IntVect(i,j,k) << " so stuck with sum_vol " << sum_vol << std::endl;
#endif
        amrex::Abort("Couldnt merge with enough cells to raise volume greater than target_volfrac");
    }
}

void
Redistribution::newlyUncoveredNbhd ( int i, int j, int k,
                                     Array4<Real const> const& apx,
                                     Array4<Real const> const& apy,
                                     Array4<Real const> const& apz,
                                     Array4<Real const> const& vfrac,
                                     Array4<Real const> const& vel_eb,
                                     Array4<int> const& itracker,
                                     Geometry const& geom,
                                     Real target_volfrac)
{
#if 0
    bool debug_print = false;
#endif

    using defaults::small_norm_diff;

    const Box domain = geom.Domain();
    const auto& is_periodic_x = geom.isPeriodic(0);
    const auto& is_periodic_y = geom.isPeriodic(1);
    const auto& is_periodic_z = geom.isPeriodic(2);

    Real apnorm, apnorm_inv;
    const Real dapx = apx(i+1,j  ,k  ) - apx(i,j,k);
    const Real dapy = apy(i  ,j+1,k  ) - apy(i,j,k);
    const Real dapz = apz(i  ,j  ,k+1) - apz(i,j,k);
    apnorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
    apnorm_inv = 1.0/apnorm;
    Real nx = dapx * apnorm_inv;
    Real ny = dapy * apnorm_inv;
    Real nz = dapz * apnorm_inv;

    bool nx_eq_ny = ( (std::abs(nx-ny) < small_norm_diff) ||
                      (std::abs(nx+ny) < small_norm_diff)  ) ? true : false;
    bool nx_eq_nz = ( (std::abs(nx-nz) < small_norm_diff) ||
                      (std::abs(nx+nz) < small_norm_diff)  ) ? true : false;
    bool ny_eq_nz = ( (std::abs(ny-nz) < small_norm_diff) ||
                      (std::abs(ny+nz) < small_norm_diff)  ) ? true : false;

    // This is a newly uncovered cell, so my vel_eb (which is time lagged here at n)
    // would be zero
    // Get an estimate from neighbors
    Real vx = 0.; // vel_eb(i,j,k,0);
    Real vy = 0.; // vel_eb(i,j,k,1);
    Real vz = 0.;

    for(int kk(-1); kk<=1; kk++) {
        for(int jj(-1); jj<=1; jj++) {
            for(int ii(-1); ii<=1; ii++) {
                if ( vfrac(i+ii,j+jj,k+kk) > 0.0 )
                {
                    if ( !(ii==0 && jj==0 && kk==0)
                         &&(is_periodic_x || (domain.smallEnd(0)<=i+ii && i+ii<=domain.bigEnd(0)))
                         &&(is_periodic_y || (domain.smallEnd(1)<=j+jj && j+jj<=domain.bigEnd(1)))
                         && Box(vel_eb).contains(IntVect(i+ii,j+jj,k+kk)) )
                    {
                        // take an "average" over all valid cells in 3x3 neighborhood
                        // excluding i,j itself. Just looking for direction, so don't worry
                        // about division
                        vx += vel_eb(i+ii,j+jj,k+kk,0);
                        vy += vel_eb(i+ii,j+jj,k+kk,1);
                        vz += vel_eb(i+ii,j+jj,k+kk,2);
//fixme
                        // if ( i==16 && j==8 ){
                        //  Print()<<"Including cell ("<<ii<<","<<jj<<"). label: "
                        //         <<label[ii+1][jj+1]<<std::endl;
                        // }
                    }
                }
            }
        }
    }

    bool vx_eq_vy = ( (std::abs(vx-vy) < small_norm_diff) ||
                      (std::abs(vx+vy) < small_norm_diff)  ) ? true : false;
    bool vx_eq_vz = ( (std::abs(vx-vz) < small_norm_diff) ||
                      (std::abs(vx+vz) < small_norm_diff)  ) ? true : false;
    bool vy_eq_vz = ( (std::abs(vy-vz) < small_norm_diff) ||
                      (std::abs(vy+vz) < small_norm_diff)  ) ? true : false;

    if ( vx_eq_vy && nx_eq_ny &&
         vx_eq_vz && nx_eq_nz &&
         vy_eq_vz && ny_eq_nz ){
        Print()<<"cell "<<IntVect(i,j,k)
               <<"eb_vel  "<< vx <<", "<<vy <<", "<<vz<<std::endl;
        Print()<<"eb_norm "<< nx <<", "<<ny <<", "<<nz<<std::endl;
        //Abort("Newly uncovered cell nbhd: Direction to merge not well resolved");
        Print() << "WARNING: Newly uncovered cell nbhd: Direction to merge not well resolved" << std::endl;
    }

    // Select first for EB motion. If that's indeterminate, choose based on EB normal
    // x-component of normal is greatest
    if ( (std::abs(vx) > std::abs(vy) && std::abs(vx) > std::abs(vz)) ||
         (std::abs(nx) > std::abs(ny) && std::abs(nx) > std::abs(nz)) )
    {
        if (nx > 0)
            itracker(i,j,k,1) = 5;
        else
            itracker(i,j,k,1) = 4;

        // y-component of normal is greatest
    } else if ( (std::abs(vy) >= std::abs(vx) && std::abs(vy) > std::abs(vz)) ||
                (std::abs(ny) >= std::abs(nx) && std::abs(ny) > std::abs(nz)) )
    {
        if (ny > 0)
            itracker(i,j,k,1) = 7;
        else

            itracker(i,j,k,1) = 2;
        // z-component of normal is greatest
    } else {
        if (nz > 0)
            itracker(i,j,k,1) = 22;
        else
            itracker(i,j,k,1) = 13;
    }

    // Override above logic if trying to reach outside a domain boundary
    // (and non-periodic)
    bool xdir_mns_ok = (is_periodic_x || (i > domain.smallEnd(0)));
    bool xdir_pls_ok = (is_periodic_x || (i < domain.bigEnd(0)  ));
    bool ydir_mns_ok = (is_periodic_y || (j > domain.smallEnd(1)));
    bool ydir_pls_ok = (is_periodic_y || (j < domain.bigEnd(1)  ));
    bool zdir_mns_ok = (is_periodic_z || (k > domain.smallEnd(2)));
    bool zdir_pls_ok = (is_periodic_z || (k < domain.bigEnd(2)  ));

    if ( (!xdir_mns_ok && (itracker(i,j,k,1) == 4)) ||
         (!xdir_pls_ok && (itracker(i,j,k,1) == 5)) )
    {
        if ( (std::abs(ny) > std::abs(nz)) )
            itracker(i,j,k,1) = (ny > 0) ? 7 : 2;
        else
            itracker(i,j,k,1) = (nz > 0) ? 22 : 13;
    }

    if ( (!ydir_mns_ok && (itracker(i,j,k,1) == 2)) ||
         (!ydir_pls_ok && (itracker(i,j,k,1) == 7)) )
    {
        if ( (std::abs(nx) > std::abs(nz)) )
            itracker(i,j,k,1) = (nx > 0) ? 5 : 4;
        else
            itracker(i,j,k,1) = (nz > 0) ? 22 : 13;
    }

    if ( (!zdir_mns_ok && (itracker(i,j,k,1) == 13)) ||
         (!zdir_pls_ok && (itracker(i,j,k,1) == 22)) )
    {
        if ( (std::abs(nx) > std::abs(ny)) )
            itracker(i,j,k,1) = (nx > 0) ? 5 : 4;
        else
            itracker(i,j,k,1) = (ny > 0) ? 7 : 2;
    }

    // (i,j,k) merges with at least one cell now
    itracker(i,j,k,0) += 1;

    // (i+ioff,j+joff,k+koff) is now the first cell in the nbhd of (i,j,k)
    auto map = getCellMap();
    int ioff = map[0][itracker(i,j,k,1)];
    int joff = map[1][itracker(i,j,k,1)];
    int koff = map[2][itracker(i,j,k,1)];

    // Sanity check
    if (vfrac(i+ioff,j+joff,k+koff) == 0.)
    {
        amrex::Print() << "Cell " << IntVect(i,j,k) << " is trying to merge with cell " << IntVect(i+ioff,j+joff,k+koff) << std::endl;
        amrex::Abort("newlyUncoveredNbhd(): Trying to merge with covered cell");
    }

    // For now, require we merge with only one other cell.
    Real sum_vol = vfrac(i,j,k) + vfrac(i+ioff,j+joff,k+koff);
    if (sum_vol < target_volfrac)
    {
        amrex::Print() << "newlyUncoveredNbhd(): Couldn't merge with enough cells to raise volume at "
                       << IntVect(i,j,k) << "to target " << target_volfrac
                       << ". Stuck with sum_vol " << sum_vol << std::endl;
        Abort("Remove abort in hydro_create_itracker_3d.cpp and recompile to continue with reduced volfrac");
    }

#if 0
    if (debug_print)
        amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vfrac(i,j,k) <<
            " trying to merge with " << IntVect(i+ioff,j+joff,k+koff) <<
            " with volfrac " << vfrac(i+ioff,j+joff,k+koff) <<
            " to get new sum_vol " <<  sum_vol << std::endl;
#endif

}
#endif
/** @} */
