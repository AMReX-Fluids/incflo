#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#include <eb_if.H>
#include <incflo.H>

using namespace amrex;

void incflo::make_eb_regular()
{
    EB2::AllRegularIF my_regular;
    auto gshop = EB2::makeShop(my_regular);
    // geom.back() is the finest AMR level; requiring max_level coarsenings
    // makes AMReX build EB data (including domain ghost cells) for every
    // AMR level, as documented for EB2::Build.
    EB2::Build(gshop, geom.back(), max_level, max_level + 100);
}
