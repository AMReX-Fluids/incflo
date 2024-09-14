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
incflo::update_vof_density (Vector<MultiFab*> const& density,Vector<MultiFab*> const& tracer)
{
  for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*tracer[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
     {
        Box const& bx = mfi.growntilebox(1);
        Array4<Real> const& density_arr = density[lev]->array(mfi);
        Array4<Real const> const& tracer_arr = tracer[lev]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {  //fixme: we use the property of the tracer 0.
           density_arr(i,j,k) = m_ro_0*(1.-tracer_arr(i,j,k,0))+m_ro_s[0]*tracer_arr(i,j,k,0);
        });
     }
	 //fixme: BCs
	density[lev]->FillBoundary(geom[lev].periodicity()); 
if(0){    	
    const auto& ba = density[lev]->boxArray();
    const auto& dm = density[lev]->DistributionMap();
    const auto& fact = density[lev]->Factory();
    // store the nodal values (the last component stores the node-centered VOF)
    MultiFab node_val(amrex::convert(ba,IntVect::TheNodeVector()),dm, 1, 0 , MFInfo(), fact);
    MultiFab center_val(ba,dm,1,0,MFInfo(), fact);	
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*density[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {        
	   Box const& nbx = surroundingNodes(mfi.tilebox());
	   Array4<Real > const& nv = node_val.array(mfi);
       Array4<Real const> const& rho   = density[lev]->const_array(mfi);	   
       ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {	
  		 
		 // calculate the node-centered VOF 
		 nv(i,j,k,0)=0.;
		 int nrho=0;		 
		 for (int detk = 0; detk > -2; --detk) 
           for (int detj = 0; detj > -2; --detj) 
             for (int deti = 0; deti > -2; --deti) {
		       Array<int,3> in{i+deti,j+detj,k+detk};
			   //averaging density to nodes
		  	   nv(i,j,k,0)+= rho(in[0],in[1],in[2],0);
		  	   nrho++; 
		       //}			 
		     }
		 nv(i,j,k,0)/= Real(nrho);		 
	   });
     }	   	   
	 average_node_to_cellcenter(center_val, 0, node_val, 0, 1);
     MultiFab::Copy(*density[lev], center_val , 0, 0, 1, 0); 
	 //fixme: BCs
	 density[lev]->FillBoundary(geom[lev].periodicity());
   } 	 
  }
	   
}


VolumeOfFluid*  incflo::get_volume_of_fluid ()
{
    if (!p_volume_of_fluid) p_volume_of_fluid = std::make_unique<VolumeOfFluid>(this);
    return p_volume_of_fluid.get();
}
