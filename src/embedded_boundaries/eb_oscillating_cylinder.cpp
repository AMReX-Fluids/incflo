#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create a simple oscillating cylinder EB.                         *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_oscillating_cylinder(Real cur_time)
{
    // Initialise cylinder parameters
    bool inside = true;
    Real radius = 0.0002;
    int direction = 0;
    Vector<Real> centervec(3);
    Vector<Real> frequency(3);
    Vector<Real> period(3);
    Vector<Real> amplitude(3);
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

    ParmParse eb("eb_flow");
    eb.getarr("frequency", frequency, 0, 3);
    eb.getarr("amplitude", amplitude, 0, 3);

    /*
       --- Adjust frequency to period ---
        The idea here is that if the frequency is too small, we just neglect it by
        setting the period to 1 and amplitude to 0.
    */
    if (std::abs(frequency[0]) > 1.e-6)
        period[0] = 1./frequency[0];
    else {
        period[0] = 1.;
        amplitude[0] = 0.;
    }

    if (std::abs(frequency[1]) > 1.e-6)
        period[1] = 1./frequency[1];
    else {
        period[1] = 1.;
        amplitude[1] = 0.;
    }

    if (std::abs(frequency[2]) > 1.e-6)
        period[2] = 1./frequency[2];
    else {
        period[2] = 1.;
        amplitude[2] = 0.;
    }

    Real PI = 3.1415926535897932384;

    Array<Real, AMREX_SPACEDIM> center = {AMREX_D_DECL(centervec[0] + amplitude[0]-amplitude[0]*cos(2*PI*frequency[0]*cur_time),
                                                       centervec[1] + amplitude[1]-amplitude[1]*cos(2*PI*frequency[1]*cur_time),
                                                       centervec[2] + amplitude[2]-amplitude[2]*cos(2*PI*frequency[2]*cur_time))};

    // FIXME:: It might be a good idea to check which direction the rotation axis is
    // Print frequency and amplitude info
#if (AMREX_SPACEDIM == 2)
    amrex::Print() << "EB Frequency: " << frequency[0] << ", "
                                       << frequency[1] << std::endl;
    amrex::Print() << "EB Amplitude: " << amplitude[0] << ", "
                                       << amplitude[1] << std::endl;
#else
    amrex::Print() << "EB Frequency: " << frequency[0] << ", "
                                       << frequency[1] << ", "
                                       << frequency[2] << std::endl;
    amrex::Print() << "EB Amplitude: " << amplitude[0] << ", "
                                       << amplitude[1] << ", "
                                       << amplitude[2] << std::endl;
#endif

    rotation = (rotation/180.)*M_PI;

    // Print info about cylinder
    amrex::Print() << " " << std::endl;
    amrex::Print() << " Internal Flow: " << inside << std::endl;
    amrex::Print() << " Radius:    " << radius << std::endl;
    amrex::Print() << " Direction: " << direction << std::endl;
    amrex::Print() << " Rotation angle(rad): " << rotation << std::endl;
    amrex::Print() << " Rotation axe: " << rotation_axe << std::endl;
#if (AMREX_SPACDEIM == 3)
    amrex::Print() << " Center:    " << center[0] << ", " << center[1] << ", " << center[2]
                   << std::endl;
#else
    amrex::Print() << " Center:    " << center[0] << ", " << center[1] << std::endl;
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
