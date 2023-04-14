#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IF_Scale.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_ParmParse.H>

#include <incflo.H>
#include <csg.hpp>

using namespace amrex;

void incflo::make_eb_csg(const std::string& geom_file, Real cur_time)
{
    bool is_internal_flow = true;
    Vector<Real> scaling_factor_vec(AMREX_SPACEDIM, 1.0);
    Vector<Real> translation_vec(AMREX_SPACEDIM, 0.0);
    Vector<Real> vels(AMREX_SPACEDIM);
    Vector<Real> omega(3);
    Vector<Real> center_of_rotation_vec(AMREX_SPACEDIM);

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
    ParmParse eb("eb_flow");
    eb.queryarr("velocity", vels, 0, AMREX_SPACEDIM);

    Array<Real,AMREX_SPACEDIM> translation;
    AMREX_D_TERM(translation[0] = translation_vec[0] + vels[0]*cur_time;,
                 translation[1] = translation_vec[1] + vels[1]*cur_time;,
                 translation[2] = translation_vec[2] + vels[2]*cur_time;);

    eb.queryarr("omega", omega, 0, 3);
    Real alpha  = omega[2]*cur_time;
    Real beta   = omega[1]*cur_time;
    Real gamma  = omega[0]*cur_time;

    Array<Real,AMREX_SPACEDIM> center_of_rotation;
    eb.queryarr("center_of_rotation", center_of_rotation_vec, 0, AMREX_SPACEDIM);
    AMREX_D_TERM(center_of_rotation[0] = center_of_rotation_vec[0];,
                 center_of_rotation[1] = center_of_rotation_vec[1];, 
                 center_of_rotation[2] = center_of_rotation_vec[2];);

    amrex::Print() << "\n Building geometry with is_internal_flow:  " << is_internal_flow << std::endl;
    auto csg_if = csg::get_csgif(geom_file, is_internal_flow);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(csg_if, "Unable to create CsgIF from geometry file");

    auto translated_csg_if = EB2::translate(EB2::scale(*csg_if, scaling_factor), translation);

    Array<Real,AMREX_SPACEDIM> to_cor;
    AMREX_D_TERM(to_cor[0] = -center_of_rotation[0];,
                 to_cor[1] = -center_of_rotation[1];, 
                 to_cor[2] = -center_of_rotation[2];);

    auto to_cor_csg_if = EB2::translate(translated_csg_if, to_cor);

    auto rotated_csg_if = EB2::rotate( 
                           EB2::rotate(
                            EB2::rotate(to_cor_csg_if, alpha, 2),
                           beta, 1), 
                          gamma, 0);

    auto final_csg_if = EB2::translate(rotated_csg_if, center_of_rotation);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(final_csg_if);

    // Build index space
    int max_level_here = 0;
    int max_coarsening_level = 100;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
