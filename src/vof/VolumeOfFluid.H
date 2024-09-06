// Header file (VolumeOfFluid class)
#ifndef INCFLO_VOF_
#define INCFLO_VOF_

class incflo;

class VolumeOfFluid
{
public:
    VolumeOfFluid (incflo* a_incflo);
    void tracer_vof_advection(amrex::Vector<amrex::MultiFab*> const& tracer,
                              AMREX_D_DECL(amrex::Vector<amrex::MultiFab const*> const& u_mac,
                                           amrex::Vector<amrex::MultiFab const*> const& v_mac,
                                           amrex::Vector<amrex::MultiFab const*> const& w_mac),
                              amrex::Real dt);
    void tracer_vof_init_fraction(int lev, amrex::MultiFab& a_tracer);
    void tracer_vof_update(int lev, amrex::MultiFab & a_tracer, amrex::Array<amrex::MultiFab, 2> & a_height);
    void curvature_calculation(int lev, amrex::MultiFab & a_tracer, amrex::Array<amrex::MultiFab, 2> & a_height,
                           amrex::MultiFab & a_kappa);
    void write_tecplot_surface(amrex::Real time, int nstep);
    void WriteTecPlotFile     (amrex::Real time, int nstep);
    void output_droplet (amrex::Real time, int nstep);
    void apply_velocity_field(amrex::Real time, int nstep);
    void velocity_face_source(int lev,amrex::Real dt, AMREX_D_DECL(amrex::MultiFab& u_mac, amrex::MultiFab& v_mac,
                                                   amrex::MultiFab& w_mac));

// normal vector of interface
   amrex::Vector<amrex::MultiFab> normal;
// the plane defined by @m.@x = @alpha
   amrex::Vector<amrex::MultiFab> alpha;
// stores the height values for calculation of curvature
// note the height is the distance between the cetroid of the current
//      and the interface in the axis direction (x-, y-, z-). Therefore,
//      each element in 'height' is a MultiFab of three components.
//      height[lev][0] stores the height values in the positive direction of the axis.
//                     it means the empty cell on the top and full cell on the bottom.
//      height[lev][1] stores the height values in the negative direction of the axis
//                     it means the empty cell on the bottom and full cell on the top.
   amrex::Vector<amrex::Array<amrex::MultiFab, 2>> height;
// Curvature
   amrex::Vector<amrex::MultiFab> kappa;
// the tag for droplets
   amrex::Vector<amrex::MultiFab> tag;

   amrex::Vector<amrex::MultiFab> force;
private:
incflo* v_incflo; //incflo object
int finest_level;
int n_tag = 0;
};



#endif
