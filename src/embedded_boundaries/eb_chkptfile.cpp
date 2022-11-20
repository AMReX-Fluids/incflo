#include <AMReX_EB2.H>

#include <incflo.H>

using namespace amrex;

void incflo::make_eb_chkptfile()
{
   // Build index space
   int max_level_here = 0;
   int max_coarsening_level = 100;
   EB2::BuildFromChkptFile("geom_chk", geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
