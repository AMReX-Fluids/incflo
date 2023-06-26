#include <hydro_redistribution.H>

using namespace amrex;

namespace Redistribution {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Array<amrex::Array<int,9>,AMREX_SPACEDIM>
getCellMap()
{
    //
    // Note that itracker has XX components and all are initialized to zero
    // We will add to the first component every time this cell is included in a merged neighborhood,
    //    either by merging or being merged
    // We identify the cells in the remaining three components with the following ordering
    //
    // ^  6 7 8
    // |  4   5
    // j  1 2 3
    //   i --->
    //
    amrex::Array<int,9> imap{0,-1,0,1,-1,1,-1,0,1};
    amrex::Array<int,9> jmap{0,-1,-1,-1,0,0,1,1,1};

    amrex::Array<amrex::Array<int,9>,AMREX_SPACEDIM> map{imap,jmap};
    return map;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Array<int,9>
getInvCellMap()
{
    //
    // Create the inverse cell map for
    //
    // ^  6 7 8
    // |  4   5
    // j  1 2 3
    //   i --->
    //
    Array<int,9> nmap{0,8,7,6,5,4,3,2,1};
    return nmap;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
normalMerging ( int i, int j, int /*k*/,
                Array4<Real const> const& apx,
                Array4<Real const> const& apy,
                Array4<Real const> const& vfrac,
                Array4<int> const& itracker,
                Geometry const& geom,
                Real target_volfrac,
                const Box& domain,
                bool is_periodic_x,
                bool is_periodic_y)
{
    int debug_verbose = 0;
    //
    // Choose based on EB normal.
    //
    auto map = Redistribution::getCellMap();
    Array<int,9> imap = map[0];
    Array<int,9> jmap = map[1];

    using defaults::small_norm_diff;

    Real apnorm, apnorm_inv, dapx, dapy;
    int k = 0;

    dapx = apx(i+1,j,k) - apx(i,j,k);
    dapy = apy(i,j+1,k) - apy(i,j,k);

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

    Real sum_vol = vfrac(i,j,k) + vfrac(i+ioff,j+joff,k);

    // if ( i==9 && j==9 ) //( debug_verbose > 0 )
    //     amrex::Print() << "Cell " << IntVect(i,j) << " with vfrac " << vfrac(i,j,k) <<
    //         " merge " << IntVect(i+ioff,j+joff) <<
    //         " with vfrac " << vfrac(i+ioff,j+joff,k) <<
    //         " to get new vfrac " <<  sum_vol << std::endl;

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

            sum_vol += vfrac(i+ioff2,j+joff2,k);

            if (debug_verbose > 0 )
            {
                std::printf("Cell (%i,%i) with volfrac %.15f trying to ALSO1 merge with (%i,%i) "
                            "with volfrac %.15f to get new sum_vol %.15f\n", 
                            i, j, vfrac(i,j,k), i+ioff2, j+joff2, vfrac(i+ioff2,j+joff2,k), sum_vol);
            }
        }
    }

    // Now we merge in the corner direction if we have already claimed two
    if (itracker(i,j,k,0) == 2)
    {
        // We already have two offsets, and we know they are in different directions
        // don't shadow
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

        sum_vol += vfrac(i+ioff,j+joff,k);

        if (debug_verbose > 0 )
        {
            std::printf("Cell (%i,%i) with volfrac %.15f trying to ALSO2 merge with (%i,%i) "
                        "with volfrac %.15f to get new sum_vol %.15f\n",
                        i, j, vfrac(i,j,k), i+ioff, j+joff, vfrac(i+ioff,j+joff,k), sum_vol);
        }
    }

    if (sum_vol < target_volfrac)
    {
        std::printf("normalMerging(): Couldn't merge with enough cells to raise volume at "
                    "(%i,%i) so stuck with sum_vol %f\n", i, j, sum_vol);
        amrex::Warning("Couldn't merge with enough cells to raise volume greater than target_volfrac");
// FIXME -- need this to get plane_right to work (has wall BC), but 0.5 is likely a better choice here
        // or using target_volume...
        if (sum_vol < Real(0.4))
        {
            amrex::Abort("Couldn't merge with enough cells to raise volume greater than 0.4");
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
newlyUncoveredNbhd ( int i, int j, int /*k*/,
                     Array4<Real const> const& apx,
                     Array4<Real const> const& apy,
                     Array4<Real const> const& vfrac,
                     Array4<Real const> const& vel_eb,
                     Array4<int> const& itracker,
                     Geometry const& geom,
                     Real target_volfrac,
                     const Box& domain,
                     bool is_periodic_x,
                     bool is_periodic_y)
{
    int debug_verbose = 1;
    //
    // Choose nbhd based on the motion of the EB
    //
    auto map = Redistribution::getCellMap();
    Array<int,9> imap = map[0];
    Array<int,9> jmap = map[1];

    using defaults::small_norm_diff;

    Real apnorm, apnorm_inv, dapx, dapy;
    int k = 0;

    dapx = apx(i+1,j,k) - apx(i,j,k);
    dapy = apy(i,j+1,k) - apy(i,j,k);

    apnorm = std::sqrt(dapx*dapx+dapy*dapy);
    apnorm_inv = 1.0/apnorm;
    Real nx = dapx * apnorm_inv;
    Real ny = dapy * apnorm_inv;

    bool nx_eq_ny = ( (std::abs(nx-ny) < small_norm_diff) ||
                      (std::abs(nx+ny) < small_norm_diff)  ) ? true : false;

    // This is a newly uncovered cell, so my vel_eb (which is time lagged here at n)
    // would be zero
    // fill with average of neighbors
    Real vx = 0.; // vel_eb(i,j,k,0);
    Real vy = 0.; // vel_eb(i,j,k,1);

    int kk=0;
    for(int jj(-1); jj<=1; jj++) {
        for(int ii(-1); ii<=1; ii++) {
            if ( vfrac(i+ii,j+jj,k) > 0.0 )
            {
                if ( !(ii==0 && jj==0)
                     &&(is_periodic_x || (domain.smallEnd(0)<=i+ii && i+ii<=domain.bigEnd(0)))
                     &&(is_periodic_y || (domain.smallEnd(1)<=j+jj && j+jj<=domain.bigEnd(1)))
                     && Box(vel_eb).contains(IntVect(i+ii,j+jj)) )
                {
                    // take an "average" over all valid cells in 3x3 neighborhood
                    // excluding i,j itself. Just looking for direction, so don't worry
                    // about division
                    vx += vel_eb(i+ii,j+jj,k+kk,0);
                    vy += vel_eb(i+ii,j+jj,k+kk,1);
//fixme
                    // if ( i==16 && j==8 ){
                    //  Print()<<"Including cell ("<<ii<<","<<jj<<"). label: "
                    //         <<label[ii+1][jj+1]<<std::endl;
                    // }
                }
            }
        }
    }

    bool vx_eq_vy = ( (std::abs(vx-vy) < small_norm_diff) ||
                      (std::abs(vx+vy) < small_norm_diff)  ) ? true : false;

    if ( vx_eq_vy && nx_eq_ny ){
        std::printf("cell (%i,%i) eb_vel %f, %f\neb_norm %f, %f\n", i, j, vx, vy, nx, ny);
        //Abort("Newly uncovered cell nbhd: Direction to merge not well resolved");
        std::printf("WARNING: Newly uncovered cell nbhd: Direction to merge not well resolved\n");
    }

    // Select first for EB motion. If that's indeterminate, choose based on EB normal
    // if (std::abs(vx) > std::abs(vy) || ( vx_eq_vy &&
    if ( std::abs(nx) > std::abs(ny) )
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

    Real sum_vol = vfrac(i,j,k) + vfrac(i+ioff,j+joff,k);

    // if ( i==10 && j==9 ) //( debug_verbose > 0 )
    //     amrex::Print() << "Cell " << IntVect(i,j) << " with vfrac " << vfrac(i,j,k) <<
    //         " merge " << IntVect(i+ioff,j+joff) <<
    //         " with vfrac " << vfrac(i+ioff,j+joff,k) <<
    //         " to get new vfrac " <<  sum_vol << std::endl;

    // For now, require we merge with only one other cell.
    if (sum_vol < target_volfrac)
    {
        std::printf("newlyUncoveredNbhd(): Couldn't merge with enough cells "
                    "to raise volume at (%i,%i) to target %f. Stuck with sum_vol %f\n", 
                    i, j, target_volfrac, sum_vol);
    }

}

//
// Central Merging
//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
centralMerging ( int i, int j,
                 Array4<Real const> const& apx,
                 Array4<Real const> const& apy,
                 Array4<Real const> const& vfrac,
                 Array4<int> const& itracker,
                 Geometry const& lev_geom,
                 Real target_volfrac,
                 const Box& domain,
                 bool is_periodic_x,
                 bool is_periodic_y)
{
    int k=0, kk=0;
    Real sum_vol = vfrac(i,j,k);

    //int label[3][3] = {{1,2,3}, {4,0,5}, {6,7,8}};
    int label[3][3] = {{1,4,6}, {2,0,7}, {3,5,8}};
    for(int jj(-1); jj<=1; jj++) {
        for(int ii(-1); ii<=1; ii++) {
            if (vfrac(i+ii,j+jj,k+kk) > 0.0)
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
                    sum_vol += vfrac(i+ii,j+jj,k+kk);

//fixme
                    if ( i==16 && j==8 ){
                        std::printf("Including cell (%i,%i). label: %f\n", ii, jj, label[ii+1][jj+1]);
                    }
                }
            }
        }
    }
// Add check on sum_vol...
}

}

/** @} */
