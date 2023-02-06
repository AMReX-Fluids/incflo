#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IF_Scale.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_ParmParse.H>

#include <incflo.H>
#include <csg.hpp>

using namespace amrex;

void incflo::make_eb_csg(const std::string& geom_file)
{
    bool is_internal_flow = true;
    Vector<Real> scaling_factor_vec(AMREX_SPACEDIM, 1.0);
    Vector<Real> translation_vec(AMREX_SPACEDIM, 0.0);

    ParmParse pp("csg");
    pp.query("internal_flow", is_internal_flow);
    if(pp.queryarr("scaling_factor", scaling_factor_vec, 0, AMREX_SPACEDIM)) {
      amrex::Print() << "WARNING: The implicit function magnitudes will not be scaled" << std::endl;
    }
    Array<Real,AMREX_SPACEDIM> scaling_factor;
    AMREX_D_TERM(scaling_factor[0] = scaling_factor_vec[0];,
                 scaling_factor[1] = scaling_factor_vec[1];, 
                 scaling_factor[2] = scaling_factor_vec[2];);

    pp.queryarr("translation", translation_vec, 0, AMREX_SPACEDIM);
    Array<Real,AMREX_SPACEDIM> translation;
    AMREX_D_TERM(translation[0] = translation_vec[0];,
                 translation[1] = translation_vec[1];,
                 translation[2] = translation_vec[2];);

    amrex::Print() << "\n Building geometry with is_internal_flow:  " << is_internal_flow << std::endl;
    auto csg_if = csg::get_csgif(geom_file, is_internal_flow);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(csg_if, "Unable to create CsgIF from geometry file");

    auto final_csg_if = EB2::translate(EB2::scale(*csg_if, scaling_factor), translation);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(final_csg_if);

    // Build index space
    int max_level_here = 0;
    int max_coarsening_level = 100;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
