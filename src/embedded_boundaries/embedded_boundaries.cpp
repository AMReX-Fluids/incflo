#include <AMReX_ParmParse.H>
#include <AMReX_EB2.H>

#include <string>
#include <algorithm>
#include <incflo.H>

using namespace amrex;

void incflo::MakeEBGeometry()
{
   /******************************************************************************
   * incflo.geometry=<string> specifies the EB geometry. <string> can be one of    *
   * box, cylinder, annulus, sphere, spherecube, twocylinders
   ******************************************************************************/

    ParmParse pp("incflo");

    std::string geom_type;
    std::string csg_file;
    pp.query("geometry", geom_type);
    pp.query("geometry_filename", csg_file);
    amrex::Print() << "incflo.geometry_filename: " << csg_file << "\n";

#ifndef CSG_EB
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( csg_file.empty(), "CSG Geometry defined in input deck but solver not built with CSG support!");
#endif

   /******************************************************************************
   *                                                                            *
   *  CONSTRUCT EB                                                              *
   *                                                                            *
   ******************************************************************************/

    if(geom_type == "cylinder")
    {
        amrex::Print() << "\n Building cylinder geometry." << "\n";
        make_eb_cylinder();
    }
    else if(geom_type == "box")
    {
        amrex::Print() << "\n Building box geometry." << "\n";
        make_eb_box();
    }
#if (AMREX_SPACEDIM == 3)
    else if(geom_type == "twocylinders")
    {
        amrex::Print() << "\n Building twocylinders geometry." << "\n";
        make_eb_twocylinders();
    }
    else if(geom_type == "spherecube")
    {
        amrex::Print() << "\n Building spherecube geometry." << "\n";
        make_eb_spherecube();
    }
    else if(geom_type == "tuscan")
    {
        amrex::Print() << "\n Building tuscan geometry." << "\n";
        make_eb_tuscan();
    }
#endif
    else if(geom_type == "annulus")
    {
        amrex::Print() << "\n Building annulus geometry." << "\n";
        make_eb_annulus();
    }
    else if(geom_type == "sphere")
    {
        amrex::Print() << "\n Building sphere geometry." << "\n";
        make_eb_sphere();
    }
    else if(geom_type == "jcap")
    {
        amrex::Print() << "\n Building JCAP geometry." << "\n";
        make_eb_cyl_tuscan();
    }
    else if(geom_type == "chkptfile")
    {
        make_eb_chkptfile();
    }
    else if(geom_type == "stl") {
        make_eb_stl();
    }
#ifdef CSG_EB
    else if(!csg_file.empty() || geom_type == "csg") {
        amrex::Print() << "\n Building geometry from .csg file:  " << csg_file << "\n";
        make_eb_csg(csg_file);
    }
#endif
    else if (geom_type.empty() || geom_type == "all_regular")
    {
        amrex::Print() << "\n No EB geometry declared in inputs => "
                       << " Will build all regular geometry." << "\n";
        make_eb_regular();
    }
    else
    {
        // An unrecognised name (a typo, a 3D-only geometry in a 2D build, or "csg"
        // in a build without CSG support) must not silently become a regular
        // geometry: the user asked for an embedded boundary and would get a plain
        // box with no diagnostic.
        std::string msg = "incflo.geometry = " + geom_type + " is not a known EB geometry";
#ifndef CSG_EB
        if (geom_type == "csg") {
            msg += " in this build (csg requires a build with CSG support)";
        }
#endif
#if (AMREX_SPACEDIM == 2)
        if (geom_type == "twocylinders" || geom_type == "spherecube" || geom_type == "tuscan") {
            msg += " in this build (it is only available in 3D)";
        }
#endif
        amrex::Abort(msg);
    }
    amrex::Print() << "Done making the EB geometry index space.\n" << "\n";

    if (m_write_geom_chk) {
        const auto& is = amrex::EB2::IndexSpace::top();
        const auto& eb_level = is.getLevel(geom.back());
        eb_level.write_to_chkpt_file("geom_chk", amrex::EB2::ExtendDomainFace(), amrex::EB2::max_grid_size);
    }
}
