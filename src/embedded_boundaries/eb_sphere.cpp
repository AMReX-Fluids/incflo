#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create a simple sphere EB.                                     *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_sphere(Real cur_time)
{
    // Initialise sphere parameters
    bool inside = true;
    Real radius = 0.0002;
    Vector<Real> centervec(3);

    // Get sphere information from inputs file.                               *
    ParmParse pp("sphere");
    ParmParse eb("eb_flow");

    pp.query("internal_flow", inside);
    pp.query("radius", radius);
    pp.getarr("center", centervec, 0, 3);

    Vector<Real> vels(AMREX_SPACEDIM);
    eb.queryarr("velocity", vels, 0, AMREX_SPACEDIM);

    // Print velocity info:
    amrex::Print() << "EB Velocity: " << vels[0] << ", " << vels[1] << std::endl;

    Array<Real, AMREX_SPACEDIM> center = {AMREX_D_DECL(centervec[0] + vels[0]*cur_time,
                                                       centervec[1] + vels[1]*cur_time,
                                                       centervec[2] + vels[2]*cur_time)};

    // Print info about sphere
    amrex::Print() << " " << std::endl;
    amrex::Print() << " Internal Flow: " << inside << std::endl;
    amrex::Print() << " Radius:    " << radius << std::endl;
#if (AMREX_SPACEDIM == 2)
    amrex::Print() << " Center:    " << center[0] << ", " << center[1] << std::endl;
#else
    amrex::Print() << " Center:    " << center[0] << ", " << center[1] << ", " << center[2]
                   << std::endl;
#endif

    // Build the sphere implicit function
    EB2::SphereIF my_sphere(radius, center, inside);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(my_sphere);

    // Build index space
    int max_level_here = 0;
    int max_coarsening_level = 100;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
