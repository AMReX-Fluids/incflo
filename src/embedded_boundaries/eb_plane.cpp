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
void incflo::make_eb_plane(Real cur_time)
{
    // Get box information from inputs file
    ParmParse pp("plane");
    ParmParse eb("eb_flow");

    if(geom[0].isAllPeriodic())
    {
        make_eb_regular();
    }
    else
    {
        /************************************************************************
         *                                                                      *
         * Define Plane geometry:                                               *
         *        -> box.normal  vector normal to the plane                     *
         *        -> box.point   vector storing a point on the plane            *
         * NOTE: walls are placed _outside_ domain for periodic directions.     *
         *                                                                      *
         ************************************************************************/

        Vector<Real> normal(AMREX_SPACEDIM, -1.0);
        pp.queryarr("normal", normal, 0, AMREX_SPACEDIM);

        Vector<Real> point(AMREX_SPACEDIM, 1.e-15);
        pp.queryarr("point", point, 0, AMREX_SPACEDIM);

        Vector<Real> vel(AMREX_SPACEDIM, 0.0);
        eb.queryarr("velocity", vel, 0, AMREX_SPACEDIM);

        // Put it correct container for EB2
        Array<Real, AMREX_SPACEDIM> loc;
        Array<Real, AMREX_SPACEDIM> norm;

        for (int i=0; i<AMREX_SPACEDIM; i++)
        {
            loc[i] = point[i] + vel[i] * cur_time;
            norm[i] = normal[i];
            Print()<<point[i]<<" "<<norm[i]<<std::endl;
        }


        EB2::PlaneIF plane_eb(loc, norm);

        // Generate GeometryShop
        auto gshop = EB2::makeShop(plane_eb);

        // Build index space
        int max_level_here = 0;
        int max_coarsening_level = 100;
        EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
    }
}
