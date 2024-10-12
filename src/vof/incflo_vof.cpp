#include <incflo.H>
#include <memory>

using namespace amrex;

void
incflo::tracer_vof_advection (Vector<MultiFab*> const& tracer, 
                             AMREX_D_DECL(Vector<MultiFab const*> const& u_mac,
                                          Vector<MultiFab const*> const& v_mac,
                                          Vector<MultiFab const*> const& w_mac))
{
    get_volume_of_fluid()->tracer_vof_advection(tracer, 
                                            AMREX_D_DECL(u_mac,v_mac,w_mac), m_dt);
}

void
incflo::update_vof_density (int lev, MultiFab & density, MultiFab & tracer)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
     {
        Box const& bx = mfi.growntilebox(1);
        Array4<Real> const& density_arr = density.array(mfi);
        Array4<Real const> const& tracer_arr = tracer.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {  //fixme: we use the property of the tracer 0.
           // Print()<<i<<" "<<j<<" "<<k<<"\n";
           //Print()<<density_arr(i,j,k)<<" "<<m_ro_0<<" "<<m_ro_s[0]<<" "<<tracer_arr(i,j,k)<<"\n";
           density_arr(i,j,k) = m_ro_0*(1.-tracer_arr(i,j,k,0))+m_ro_s[0]*tracer_arr(i,j,k,0);           
        });
     }
     //fixme: BCs
    density.FillBoundary(geom[lev].periodicity()); 

       
}


VolumeOfFluid*  incflo::get_volume_of_fluid ()
{
    if (!p_volume_of_fluid) p_volume_of_fluid = std::make_unique<VolumeOfFluid>(this);
    return p_volume_of_fluid.get();
}
