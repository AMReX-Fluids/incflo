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
void incflo::make_eb_hexahedron(Real cur_time)
{
    // Get box information from inputs file
    ParmParse pp("hexahedron");
    ParmParse eb("eb_flow");

    /************************************************************************
     *                                                                      *
     * Define Box geometry:                                                 *
     *        -> box.{Lo,Hi} vector storing box lo/hi                       *
     *        -> box.offset  vector storing box offset                      *
     *                                                                      *
     ************************************************************************/

    Vector<Real> boxLo(AMREX_SPACEDIM), boxHi(AMREX_SPACEDIM);
    Real offset = 1.0e-15;
    Real angle = 0.75;
    int dir = 0;
    int dir2 = 1;

    for(int i = 0; i < AMREX_SPACEDIM; i++)
    {
        boxLo[i] = geom[0].ProbLo(i);
        boxHi[i] = geom[0].ProbHi(i);
    }

    pp.queryarr("lo", boxLo, 0, AMREX_SPACEDIM);
    pp.queryarr("hi", boxHi, 0, AMREX_SPACEDIM);
    pp.query("angle", angle);
    pp.query("dir", dir);

    Vector<Real> vel(AMREX_SPACEDIM);
    eb.queryarr("velocity", vel, 0, AMREX_SPACEDIM);

    pp.query("offset", offset);


    Real xlo = boxLo[0] + offset;
    Real xhi = boxHi[0] - offset;

    Real ylo = boxLo[1] + offset;
    Real yhi = boxHi[1] - offset;

    xlo = xlo + vel[0] * cur_time;
    xhi = xhi + vel[0] * cur_time;

    ylo = ylo + vel[1] * cur_time;
    yhi = yhi + vel[1] * cur_time;

#if AMREX_SPACEDIM > 2 
    Real zlo = boxLo[2] + offset;
    Real zhi = boxHi[2] - offset;

    zlo = zlo + vel[2] * cur_time;
    zhi = zhi + vel[2] * cur_time;
#endif

    bool inside = false;
    RealArray lo {AMREX_D_DECL(xlo, ylo, zlo)};
    RealArray hi {AMREX_D_DECL(xhi, yhi, zhi)};

//    auto polyhedron = EB2::rotate(EB2::rotate(EB2::BoxIF(lo, hi, inside), angle, dir), angle, dir2);
    auto polyhedron = EB2::rotate(EB2::BoxIF(lo, hi, inside), angle, dir);
    auto gshop = EB2::makeShop(polyhedron);

    // Build index space
    int max_level_here = 0;
    int max_coarsening_level = 100;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}

