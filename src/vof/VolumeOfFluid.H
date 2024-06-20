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
    void tracer_vof_init_fraction(amrex::Vector<amrex::MultiFab*> const& tracer);		               
    void write_tecplot_surface(amrex::Real time, int nstep);
    void WriteTecPlotFile     (amrex::Real time, int nstep);		               	            	

// normal vector of interface 
   amrex::Vector<amrex::MultiFab> normal; 
// the plane defined by @m.@x = @alpha
   amrex::Vector<amrex::MultiFab> alpha;
private:
incflo* v_incflo; //incflo object
int finest_level;
#if AMREX_SPACEDIM == 2
# define F(x,y,z) f[x][y]
#else
# define F(x,y,z) f[x][y][z]
#endif

void tracer_vof_update(amrex::Vector<amrex::MultiFab*> const& tracer);
};

#endif
