#include <AMReX_EB2.H>

#include <incflo.H>

using namespace amrex;

void incflo::make_eb_chkptfile()
{
   // Build index space
   // geom.back() is the finest AMR level; requiring max_level coarsenings
   // makes AMReX build EB data (including domain ghost cells) for every
   // AMR level, as documented for EB2::Build.
   int max_level_here = max_level;
   int max_coarsening_level = 100;
   EB2::BuildFromChkptFile("geom_chk", geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
