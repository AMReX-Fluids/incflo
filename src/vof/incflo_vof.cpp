#include <incflo.H>
#include <memory>

using namespace amrex;

void
incflo::tracer_vof_advection(Vector<MultiFab*> const& tracer, 
                             AMREX_D_DECL(Vector<MultiFab const*> const& u_mac,
                                          Vector<MultiFab const*> const& v_mac,
                                          Vector<MultiFab const*> const& w_mac))
{
    get_volume_of_fluid()->tracer_vof_advection(tracer, 
	                                        AMREX_D_DECL(u_mac,v_mac,w_mac), m_dt);
}



VolumeOfFluid*
incflo::get_volume_of_fluid ()
{
    if (!p_volume_of_fluid) p_volume_of_fluid = std::make_unique<VolumeOfFluid>(this);
    return p_volume_of_fluid.get();
}
