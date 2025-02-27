// Header file (VolumeOfFluid class)
#ifndef INCFLO_VOF_
#define INCFLO_VOF_

#include <AMReX_Interpolater.H>

class incflo;

class VolumeOfFluid
{
public:
    friend incflo;
    VolumeOfFluid (incflo* a_incflo);
    void tracer_vof_advection(amrex::Vector<amrex::MultiFab*> const& tracer,
                              AMREX_D_DECL(amrex::Vector<amrex::MultiFab const*> const& u_mac,
                                           amrex::Vector<amrex::MultiFab const*> const& v_mac,
                                           amrex::Vector<amrex::MultiFab const*> const& w_mac),
                              amrex::Real dt);
    void tracer_vof_init_fraction(int lev, amrex::MultiFab& a_tracer);
    void tracer_vof_update(int lev, amrex::MultiFab const & a_tracer, amrex::Array<amrex::MultiFab, 2> & a_height);
    void curvature_calculation(int lev, amrex::MultiFab const & a_tracer, amrex::Array<amrex::MultiFab, 2> & a_height,
                           amrex::MultiFab & a_kappa);
    void write_tecplot_surface(int max_level, amrex::Real time, int nstep);
    void WriteTecPlotFile     (int max_level, amrex::Real time, int nstep);
    void output_droplet (amrex::Real time, int nstep);
     int domain_tag_droplets (int finest_level, amrex::Vector<amrex::BoxArray > const &grids, amrex::Vector<amrex::Geometry > const& geom,
                              amrex::Vector<amrex::MultiFab*> const& vof,amrex::Vector<amrex::MultiFab*> const& tag);
    void apply_velocity_field(amrex::Real time, int nstep);
    void velocity_face_source(int lev,amrex::Real dt, AMREX_D_DECL(amrex::MultiFab& u_mac, amrex::MultiFab& v_mac,
                                                                   amrex::MultiFab& w_mac),
                              AMREX_D_DECL(amrex::MultiFab* sfu_mac, amrex::MultiFab* sfv_mac, amrex::MultiFab* sfw_mac));
    void variable_filtered (int lev, amrex::MultiFab& variable);
    void curvature_average_down (amrex::MultiFab const & kappa_fine, amrex::MultiFab & kappa_crse, amrex::IntVect const & ratio);


/**
*
* VOF interpolation on cell centered data.
*/
class VOFCoarseToFine : public amrex::Interpolater
{
public:
    /**
    * Returns coarsened box given fine box and refinement ratio.
    */
    amrex::Box CoarseBox (const amrex::Box& fine, int ratio) override;

    /**
    * Returns coarsened box given fine box and refinement ratio.
    */
    amrex::Box CoarseBox (const amrex::Box& fine, const amrex::IntVect& ratio) override;

    /**
    * Coarse to fine interpolation in space.
    */
    void interp (const amrex::FArrayBox& crse, int crse_comp,
                       amrex::FArrayBox& fine, int fine_comp, int ncomp,
                 const amrex::Box& fine_region, const amrex::IntVect&  ratio,
                 const amrex::Geometry&  crse_geom, const amrex::Geometry&  fine_geom,
                 amrex::Vector<amrex::BCRec> const& bcr,
                 int              actual_comp,int actual_state,
                 amrex::RunOn runon) override;
    int lev;
};


private:
   incflo* v_incflo; //incflo object
//   int finest_level;
   int n_tag = 0;
   int output_drop_frequence = 1;

   struct LevelData {
        LevelData () = default;
        LevelData (amrex::BoxArray const& ba,
                   amrex::DistributionMapping const& dm,
                   amrex::FabFactory<amrex::FArrayBox> const& fact,
                   incflo* a_incflo);
        // cell-centered multifabs
// normal vector of interface "m"
        amrex::MultiFab normal;
// the plane defined by @m.@x = @alpha
        amrex::MultiFab alpha;
// Curvature
        amrex::MultiFab kappa;
// the tag for droplets
        amrex::MultiFab tag;
// the tag for surface force
        amrex::MultiFab force;
// stores the height values for calculation of curvature
// note the height is the distance between the cetroid of the current
//      and the interface in the axis direction (x-, y-, z-). Therefore,
//      each element in 'height' is a MultiFab of three components.
//      height[lev][0] stores the height values in the positive direction of the axis.
//                     it means the empty cell on the top and full cell on the bottom.
//      height[lev][1] stores the height values in the negative direction of the axis
//                     it means the empty cell on the bottom and full cell on the top.
        amrex::Array<amrex::MultiFab, 2> height;
   };
   amrex::Vector<std::unique_ptr<LevelData> > m_leveldata;

   amrex::Vector<amrex::MultiFab*>       get_vector_ptr(std::function<amrex::MultiFab&(LevelData&)> member_selector) noexcept;
   amrex::Vector<amrex::MultiFab const*> get_vector_const_ptr(std::function<amrex::MultiFab&(LevelData&)> member_selector) noexcept;
   VOFCoarseToFine    vof_interp;
};



#endif
