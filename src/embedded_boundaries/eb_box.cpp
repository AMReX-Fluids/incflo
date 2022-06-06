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
void incflo::make_eb_box(Real cur_time)
{
    // Get box information from inputs file
    ParmParse pp("box");
    ParmParse eb("eb_flow");

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

        for(int i = 0; i < AMREX_SPACEDIM; i++)
        {
            boxLo[i] = geom[0].ProbLo(i);
            boxHi[i] = geom[0].ProbHi(i);
        }

        pp.queryarr("Lo", boxLo, 0, AMREX_SPACEDIM);
        pp.queryarr("Hi", boxHi, 0, AMREX_SPACEDIM);  
        
        Real velx = 0.0;
        eb.query("vel_mag", velx);

        pp.query("offset", offset);

        amrex::Print() << "velx * cur_time: " << velx *cur_time << std::endl;

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

#if (AMREX_SPACEDIM == 2)

        xlo = xlo + velx * cur_time;
        xhi = xhi + velx * cur_time; 
        
        Array<Real, 2> point_lox{xlo, 0.0};
        Array<Real, 2> normal_lox{-1.0, 0.0};
        Array<Real, 2> point_hix{xhi, 0.0};
        Array<Real, 2> normal_hix{1.0, 0.0};

        Array<Real, 2> point_loy{0.0, ylo};
        Array<Real, 2> normal_loy{0.0, -1.0};
        Array<Real, 2> point_hiy{0.0, yhi};
        Array<Real, 2> normal_hiy{0.0, 1.0};

        EB2::PlaneIF plane_lox(point_lox, normal_lox);
        EB2::PlaneIF plane_hix(point_hix, normal_hix);

        EB2::PlaneIF plane_loy(point_loy, normal_loy);
        EB2::PlaneIF plane_hiy(point_hiy, normal_hiy);

        // Generate GeometryShop
        auto gshop = EB2::makeShop(EB2::makeUnion(plane_lox, plane_hix,
                                                  plane_loy, plane_hiy));
#else
        Real zlo = boxLo[2] + offset;
        Real zhi = boxHi[2] - offset;

        // This ensures that the walls won't even touch the ghost cells. By
        // putting them one domain width away
        if(geom[0].isPeriodic(2))
        {
            zlo = 2.0 * geom[0].ProbLo(2) - geom[0].ProbHi(2);
            zhi = 2.0 * geom[0].ProbHi(2) - geom[0].ProbLo(2);
        }

        Array<Real, 3> point_lox{xlo, 0.0, 0.0};
        Array<Real, 3> normal_lox{-1.0, 0.0, 0.0};
        Array<Real, 3> point_hix{xhi, 0.0, 0.0};
        Array<Real, 3> normal_hix{1.0, 0.0, 0.0};

        Array<Real, 3> point_loy{0.0, ylo, 0.0};
        Array<Real, 3> normal_loy{0.0, -1.0, 0.0};
        Array<Real, 3> point_hiy{0.0, yhi, 0.0};
        Array<Real, 3> normal_hiy{0.0, 1.0, 0.0};

        Array<Real, 3> point_loz{0.0, 0.0, zlo};
        Array<Real, 3> normal_loz{0.0, 0.0, -1.0};
        Array<Real, 3> point_hiz{0.0, 0.0, zhi};
        Array<Real, 3> normal_hiz{0.0, 0.0, 1.0};

        EB2::PlaneIF plane_lox(point_lox, normal_lox);
        EB2::PlaneIF plane_hix(point_hix, normal_hix);

        EB2::PlaneIF plane_loy(point_loy, normal_loy);
        EB2::PlaneIF plane_hiy(point_hiy, normal_hiy);

        EB2::PlaneIF plane_loz(point_loz, normal_loz);
        EB2::PlaneIF plane_hiz(point_hiz, normal_hiz);

        // Generate GeometryShop
        auto gshop = EB2::makeShop(EB2::makeUnion(plane_lox, plane_hix,
                                                  plane_loy, plane_hiy,
                                                  plane_loz, plane_hiz));
#endif

        // Build index space
        int max_level_here = 0;
        int max_coarsening_level = 100;
        EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
    }
}
