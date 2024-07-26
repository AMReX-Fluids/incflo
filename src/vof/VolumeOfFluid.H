//struct TracerVOF {
//  GfsVariableTracer parent;
//
//  GfsVariable * m[FTT_DIMENSION], * alpha;
//};
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
    void write_tecplot_surface(amrex::Real time, int nstep);
    void WriteTecPlotFile     (amrex::Real time, int nstep);
    void output_droplet (amrex::Real time, int nstep);
    void apply_velocity_field(amrex::Real time, int nstep);

// normal vector of interface
   amrex::Vector<amrex::MultiFab> normal;
// the plane defined by @m.@x = @alpha
   amrex::Vector<amrex::MultiFab> alpha;
// the heights used for calculation of curvature
   amrex::Vector<amrex::Array<amrex::MultiFab, 2>> height;
// the tag for droplets
   amrex::Vector<amrex::MultiFab> tag;
private:
incflo* v_incflo; //incflo object
int finest_level;
int n_tag = 0;
};



#endif
