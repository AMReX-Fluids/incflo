/*
 * \file hydro_create_itracker_2d.cpp
 * \addtogroup Redistribution
 * @{
 *
 */

#include <hydro_redistribution.H>

using namespace amrex;

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
    int debug_verbose = 1;

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

    amrex::ParallelFor(bx_per_g4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We check for cut-cells in the new geometry
        if ( (vfrac_new(i,j,k) > 0.0 && vfrac_new(i,j,k) < 1.0) && vfrac_old(i,j,k) > 0.0)
        {
            normalMerging(i, j, k,
			  AMREX_D_DECL(apx_new, apy_new, apz_new),
			  vfrac_new, itracker,
                          lev_geom, target_volfrac);
        }
        else if ( (vfrac_new(i,j,k) > 0.0 && vfrac_new(i,j,k) < 1.0) && vfrac_old(i,j,k) == 0.0)
        {
            // For now, require that newly uncovered cells only have one other cell in it's nbhd
            // FIXME, unsure of target_volfrac here...
            newlyUncoveredNbhd(i, j, k,
			       AMREX_D_DECL(apx_new, apy_new, apz_new),
			       vfrac_new, vel_eb, itracker,
                               lev_geom, 0.5);
        }
        else if ( vfrac_old(i,j,k) > 0.0 && vfrac_new(i,j,k) == 0.0)
        {
            // Create a nbhd for cells that become covered...
            // vfrac is only for checking volume of nbhd
            // Probably don't need target_volfrac to match with general case,
            // only need to put this in one cell???
            normalMerging(i, j, k,
			  AMREX_D_DECL(apx_new, apy_new, apz_new),
			  vfrac_new, itracker,
                          lev_geom, target_volfrac);
        }
    });


#if 0
    amrex::Print() << "\nInitial Cell Merging" << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
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
    });
    amrex::Print() << std::endl;
#endif

    // FIXME - Need to do some sort of check for if they've already been included though
    // Check uncovered and covered cells, make sure the neighbors also include them.
    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Newly covered Cells
        if (vfrac_new(i,j,k) == 0. && vfrac_old(i,j,k) > 0.0)
        {
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
		int ioff = map[0][itracker(i,j,k,i_nbor)];
		int joff = map[1][itracker(i,j,k,i_nbor)];
		int koff = (AMREX_SPACEDIM < 3) ? 0 : map[2][itracker(i,j,k,i_nbor)];

                if ( Box(itracker).contains(Dim3{i+ioff,j+joff,k+koff}) )
                {
                    // amrex::Print() << "Cell  " << Dim3{i,j,k} << " is covered and merged with neighbor at " << Dim3{i+ioff,j+joff,k+koff} << std::endl;
                    itracker(i+ioff,j+joff,k+koff,0) += 1;
		    int comp = itracker(i+ioff,j+joff,k+koff,0);
                    itracker(i+ioff,j+joff,k+koff,comp) = nmap[itracker(i,j,k,i_nbor)];
                }
            }
        }

        // Newly uncovered
        if (vfrac_new(i,j,k) > 0. && vfrac_new(i,j,k) < 1. && vfrac_old(i,j,k) == 0.0 )
        {
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
		int ioff = map[0][itracker(i,j,k,i_nbor)];
		int joff = map[1][itracker(i,j,k,i_nbor)];
		int koff = (AMREX_SPACEDIM < 3) ? 0 : map[2][itracker(i,j,k,i_nbor)];

                if ( Box(itracker).contains(Dim3{i+ioff,j+joff,k+koff}) )
                {
                    // amrex::Print() << "Cell  " << Dim3{i,j,k} << " is newly uncovered and merged with neighbor at " << Dim3{i+ioff,j+joff,k+koff} << std::endl;
                    itracker(i+ioff,j+joff,k+koff,0) += 1;
		    int comp = itracker(i+ioff,j+joff,k+koff,0);
                    itracker(i+ioff,j+joff,k+koff,comp) = nmap[itracker(i,j,k,i_nbor)];
                }
            }
        }
    });

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

#if 1
    amrex::Print() << "Post Update to Cell Merging" << std::endl;

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
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
    });
    amrex::Print() << std::endl;
#endif

}
/** @} */
