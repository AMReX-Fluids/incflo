#include <AMReX_ParmParse.H>
#include <AMReX_EB2.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

void incflo::MakeEBGeometry(Real cur_time)
{
   /******************************************************************************
   * incflo.geometry=<string> specifies the EB geometry. <string> can be one of    *
   * box, cylinder, annulus, sphere, spherecube, twocylinders
   ******************************************************************************/

    ParmParse pp("incflo");

    std::string geom_type;
    pp.query("geometry", geom_type);

   /******************************************************************************
   *                                                                            *
   *  CONSTRUCT EB                                                              *
   *                                                                            *
   ******************************************************************************/

    if(geom_type == "cylinder")
    {
    amrex::Print() << "\n Building cylinder geometry." << std::endl;
        make_eb_cylinder();
    }
    else if(geom_type == "box")
    {
        amrex::Print() << "\n Building box geometry." << std::endl;
        make_eb_box(cur_time);
    }
#if (AMREX_SPACEDIM == 3)
    else if(geom_type == "twocylinders")
    {
    amrex::Print() << "\n Building twocylinders geometry." << std::endl;
        make_eb_twocylinders();
    }
    else if(geom_type == "spherecube")
    {
    amrex::Print() << "\n Building spherecube geometry." << std::endl;
        make_eb_spherecube();
    }
    else if(geom_type == "tuscan")
    {
    amrex::Print() << "\n Building tuscan geometry." << std::endl;
        make_eb_tuscan();
    }
#endif
    else if(geom_type == "annulus")
    {
    amrex::Print() << "\n Building annulus geometry." << std::endl;
        make_eb_annulus();
    }
    else if(geom_type == "sphere")
    {
    amrex::Print() << "\n Building sphere geometry at time " << cur_time << std::endl;
        make_eb_sphere(cur_time);
    }
    else if(geom_type == "stirrer")
    {
    amrex::Print() << "\n Building stirrer geometry at time " << cur_time << std::endl;
        make_eb_stirrer(cur_time);
    }
    else if(geom_type == "jcap")
    {
    amrex::Print() << "\n Building JCAP geometry." << std::endl;
        make_eb_cyl_tuscan();
    }
    else if(geom_type == "chkptfile")
    {
       make_eb_chkptfile();
    }
    else
    {
    amrex::Print() << "\n No EB geometry declared in inputs => "
                   << " Will build all regular geometry." << std::endl;
        make_eb_regular();
    }
    amrex::Print() << "Done making the geometry ebfactory.\n" << std::endl;

    if (m_write_geom_chk) {
       const auto& is = amrex::EB2::IndexSpace::top();
       const auto& eb_level = is.getLevel(geom.back());
       eb_level.write_to_chkpt_file("geom_chk", amrex::EB2::ExtendDomainFace(), amrex::EB2::max_grid_size);
    }
}
