#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

/****************************************************************************
 * Function to create a simple rectangular box with EB walls.               *
 *                                                                          *
 ****************************************************************************/
void incflo::make_eb_box()
{
    // Get box information from inputs file
    ParmParse pp("box");

    if(geom[0].isAllPeriodic())
    {
        make_eb_regular();
    }
    else
    {
        /************************************************************************
         *                                                                      *
         * Define Box geometry:                                                 *
         *        -> box.{Lo,Hi} vector storing box lo/hi                       *
         *        -> box.offset  vector storing box offset                      *
         * NOTE: walls are placed _outside_ domain for periodic directions.     *
         *                                                                      *
         ************************************************************************/

        Vector<Real> boxLo(AMREX_SPACEDIM), boxHi(AMREX_SPACEDIM);
        Real offset = 1.0e-15;
        bool inside = true;
        
        for(int i = 0; i < AMREX_SPACEDIM; i++)
        {
            boxLo[i] = geom[0].ProbLo(i);
            boxHi[i] = geom[0].ProbHi(i);
        }

        pp.queryarr("Lo", boxLo, 0, AMREX_SPACEDIM);
        pp.queryarr("Hi", boxHi, 0, AMREX_SPACEDIM);

        pp.query("offset", offset);
        pp.query("internal_flow", inside);

        Real xlo = boxLo[0] + offset;
        Real xhi = boxHi[0] - offset;

        // This ensures that the walls won't even touch the ghost cells. By
        // putting them one domain width away
        if(geom[0].isPeriodic(0))
        {
            xlo = 2.0 * geom[0].ProbLo(0) - geom[0].ProbHi(0);
            xhi = 2.0 * geom[0].ProbHi(0) - geom[0].ProbLo(0);
        }

        Real ylo = boxLo[1] + offset;
        Real yhi = boxHi[1] - offset;

        // This ensures that the walls won't even touch the ghost cells. By
        // putting them one domain width away
        if(geom[0].isPeriodic(1))
        {
            ylo = 2.0 * geom[0].ProbLo(1) - geom[0].ProbHi(1);
            yhi = 2.0 * geom[0].ProbHi(1) - geom[0].ProbLo(1);
        }

#if (AMREX_SPACEDIM > 2)
        Real zlo = boxLo[2] + offset;
        Real zhi = boxHi[2] - offset;

        // This ensures that the walls won't even touch the ghost cells. By
        // putting them one domain width away
        if(geom[0].isPeriodic(2))
        {
            zlo = 2.0 * geom[0].ProbLo(2) - geom[0].ProbHi(2);
            zhi = 2.0 * geom[0].ProbHi(2) - geom[0].ProbLo(2);
        }
#endif

        RealArray lo {AMREX_D_DECL(xlo, ylo, zlo)};
        RealArray hi {AMREX_D_DECL(xhi, yhi, zhi)};
        
        EB2::BoxIF my_box(lo, hi, inside);

        // Generate GeometryShop
        auto gshop = EB2::makeShop(my_box);
        
        // Build index space
        int max_level_here = 0;
        int max_coarsening_level = 100;
        EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
    }
}
