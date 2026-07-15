#include <AMReX_EB2.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EB2_IndexSpace_STL.H>

#include <incflo.H>

using namespace amrex;

void incflo::make_eb_stl ()
{
    bool is_internal_flow = true;
    Real scaling_factor = 1.0;
    Vector<Real> translation_vec(3, 0.0);
    bool use_bvh = true;

    ParmParse pp("stl");

    std::string stl_file;
    pp.query("geometry_filename", stl_file);

    if (stl_file.empty()) {
        Abort("\nMissing or invalid input: stl.geometry_filename = ...");
    } else {
        Print() << "STL geometry file: " << stl_file << '\n';
    }

    pp.query("internal_flow", is_internal_flow); // This will flip the normal
    pp.query("scaling_factor", scaling_factor);
    pp.queryarr("translation", translation_vec); // This is the stl center
    pp.query("use_bvh", use_bvh);

    if (is_internal_flow) { Print() << "\n Building geometry for internal flow\n"; }
    else { Print() << "\n Building geometry for external flow\n"; }

    /************************************************************************
     *                                                                      *
     * Build EB levels from STL                                             *
     *                                                                      *
     ***********************************************************************/
    EB2::IndexSpace::push(
        std::make_unique<EB2::IndexSpaceSTL>
        (stl_file, scaling_factor,
         Array<Real,3>{translation_vec[0], translation_vec[1], translation_vec[2]},
         is_internal_flow, geom[max_level], max_level, 100, /*ngrow*/4,
         /* build coarse by coarsening*/ true, EB2::ExtendDomainFace(),
         EB2::NumCoarsenOpt(), use_bvh, /*support mvmc*/false));

}
