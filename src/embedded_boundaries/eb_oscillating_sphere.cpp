#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create a simple oscillating sphere EB.                           *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_oscillating_sphere(Real cur_time)
{
    amrex::Print() << " " << std::endl;
    amrex::Print() << " Making oscillating sphere geometry: " << std::endl;

    // Initialize sphere parameters
    bool inside = true;
    Real radius = 0.0002;
    Vector<Real> centervec(3);
    Vector<Real> frequency(3);
    Vector<Real> period(3);
    Vector<Real> amplitude(3);

    // Get sphere information from inputs file.
    ParmParse pp("sphere");

    pp.query("internal_flow", inside);
    pp.query("radius", radius);
    pp.getarr("center", centervec, 0, 3);

    // Get frequency and amplitude.
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

    Array<Real, AMREX_SPACEDIM> center = {AMREX_D_DECL(centervec[0] + amplitude[0]*sin(period[0]*cur_time),
                                                       centervec[1] + amplitude[1]*sin(period[1]*cur_time),
                                                       centervec[2] + amplitude[2]*sin(period[2]*cur_time))};

    // Print frequency and amplitude info
#if (AMREX_SPACEDIM == 2)
    amrex::Print() << "EB Frequency: " << frequency[0] << ", " << frequency[1] << std::endl;
    amrex::Print() << "EB Amplitude: " << amplitude[0] << ", " << amplitude[1] << std::endl;
#else
    amrex::Print() << "EB Frequency: " << frequency[0] << ", " << frequency[1] << ", " << frequency[2] << std::endl;
    amrex::Print() << "EB Amplitude: " << amplitude[0] << ", " << amplitude[1] << ", " << amplitude[2] << std::endl;
#endif

   // Print info about the sphere
   amrex::Print() << " " << std::endl;
   amrex::Print() << " Internal Flow: " << inside << std::endl;
   amrex::Print() << " Radius: " << radius << std::endl;
#if (AMREX_SPACEDIM == 2)
   amrex::Print() << " Center: " << center[0] << ", " << center[1] << std::endl;
#else
   amrex::Print() << " Center: " << center[0] << ", " << center[1] << ", " << center[2] << std::endl;
#endif

    // Build the sphere implicit function
    EB2::SphereIF my_sphere(radius, center, inside);

    //Generate GeometryShop
    auto gshop = EB2::makeShop(my_sphere);

    // Build index space
    int max_level_here = 0;
    int max_coarsening_level = 100;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
