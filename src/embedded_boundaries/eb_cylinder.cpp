#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create a simple cylinder EB.                                     *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_cylinder()
{
    // Initialise cylinder parameters
    bool inside = true;
    Real radius = 0.0002;
    int direction = 0;
    Vector<Real> centervec(3);
    Real rotation  = 0;
    int rotation_axe  = 0;

    // Get cylinder information from inputs file.                               *
    ParmParse pp("cylinder");

    pp.query("internal_flow", inside);
    pp.query("radius", radius);
    pp.query("direction", direction);
    pp.query("rotation",   rotation);
    pp.query("rotation_axe",   rotation_axe);
    pp.getarr("center", centervec, 0, 3);
    Array<Real, AMREX_SPACEDIM> center = {AMREX_D_DECL(centervec[0], centervec[1], centervec[2])};

    rotation = (rotation/Real(180))*Real(M_PI);

    // Print info about cylinder
    amrex::Print() << " " << "\n";
    amrex::Print() << " Internal Flow: " << inside << "\n";
    amrex::Print() << " Radius:    " << radius << "\n";
    amrex::Print() << " Direction: " << direction << "\n";
    amrex::Print() << " Rotation angle(rad): " << rotation << "\n";
    amrex::Print() << " Rotation axe: " << rotation_axe << "\n";
#if (AMREX_SPACEDIM == 3)
    amrex::Print() << " Center:    " << center[0] << ", " << center[1] << ", " << center[2]
                   << "\n";
#else
    amrex::Print() << " Center:    " << center[0] << ", " << center[1] << "\n";
#endif

    // Build the Cylinder implficit function representing the curved walls
    EB2::CylinderIF my_cyl(radius, direction, center, inside);

    auto my_cyl_rot = EB2::rotate(my_cyl, rotation, rotation_axe);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(my_cyl_rot);

    // Build index space
    int max_level_here = 0;
    int max_coarsening_level = 100;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
