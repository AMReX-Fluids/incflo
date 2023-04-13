#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create a simple sphereCUBE EB.                                     *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_stirrer(Real time)
{
    amrex::Print() << " " << std::endl;
    amrex::Print() << " Making stirrer geometry: " << std::endl;

    Real domain_length =  geom[0].ProbHi(0) - geom[0].ProbLo(0);

    bool inside = false;

    Real width  = 0.1*domain_length;
    Real length = 0.4*domain_length;
    Real radius = 0.5*width;

    amrex::Print() << "Domain  length " << domain_length << std::endl;
    amrex::Print() << "Stirrer length " <<        length << std::endl;
    amrex::Print() << "Stirrer width  " <<         width << std::endl;
    amrex::Print() << "Sphere  radius " <<        radius << std::endl;

    AMREX_D_TERM(Real domain_center_x = 0.5 * (geom[0].ProbLo(0) + geom[0].ProbHi(0));,
                 Real domain_center_y = 0.5 * (geom[0].ProbLo(1) + geom[0].ProbHi(1));,
                 Real domain_center_z = 0.5 * (geom[0].ProbLo(2) + geom[0].ProbHi(2)););

    Array<Real,AMREX_SPACEDIM> center1 = {AMREX_D_DECL(domain_center_x,domain_center_y-0.5*length,domain_center_z)};
    Array<Real,AMREX_SPACEDIM> center2 = {AMREX_D_DECL(domain_center_x,domain_center_y+0.5*length,domain_center_z)};

    // Build the sphere implicit functions
    EB2::SphereIF sphere1(radius, center1, inside);
    EB2::SphereIF sphere2(radius, center2, inside);
    auto spheres = EB2::makeUnion(sphere1, sphere2);

    Array<Real,AMREX_SPACEDIM> rectangle_lo = {AMREX_D_DECL(0.5-radius,0.5-0.5*length,0.5)};
    Array<Real,AMREX_SPACEDIM> rectangle_hi = {AMREX_D_DECL(0.5+radius,0.5+0.5*length,0.5)};
    EB2::BoxIF rectangle(rectangle_lo, rectangle_hi, inside);

    auto stirrer = EB2::makeUnion(spheres, rectangle);

    // Rotate the stirrer
    ParmParse eb("eb_flow");
    Vector<Real> omega(3);
    eb.queryarr("omega", omega, 0, 3);

    Real alpha  = omega[2]*time;
    Real beta   = omega[1]*time;
    Real gamma  = omega[0]*time;

    auto stirrer_at_center = EB2::translate(stirrer, {
                                AMREX_D_DECL(-domain_center_x,-domain_center_y,-domain_center_z)});

    auto rotated_stirrer_at_center = EB2::rotate(
                                        EB2::rotate(
                                           EB2::rotate(
                                              stirrer_at_center, alpha, 2),
                                        beta, 1), gamma, 0);

    auto stirrer_final = EB2::translate(rotated_stirrer_at_center, {
                                AMREX_D_DECL(domain_center_x,domain_center_y,domain_center_z)});

    // Generate GeometryShop
    auto gshop = EB2::makeShop(stirrer_final);

    // Build index space
    int max_level_here = 0;
    int max_coarsening_level = 100;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
