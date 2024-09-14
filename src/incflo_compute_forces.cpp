#include <incflo.H>

using namespace amrex;

void incflo::compute_tra_forces (Vector<MultiFab*> const& tra_forces,
                                 Vector<MultiFab const*> const& density)
{
    // NOTE: this routine must return the force term for the update of (rho s), NOT just s.
    if (m_advect_tracer) {

        auto const* iconserv = get_tracer_iconserv_device_ptr();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (MFIter mfi(*tra_forces[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real>       const& tra_f = tra_forces[lev]->array(mfi);
                Array4<Real const> const& rho   =    density[lev]->const_array(mfi);

                ParallelFor(bx, m_ntrac,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    // For now we don't have any external forces on the scalars
                    tra_f(i,j,k,n) = 0.0;

                    if (iconserv[n]){
                        // Return the force term for the update of (rho s), NOT just s.
                        tra_f(i,j,k,n) *= rho(i,j,k);
                    }
                });
            }
        }
    }
}

void incflo::compute_vel_forces (Vector<MultiFab*> const& vel_forces,
                                 Vector<MultiFab const*> const& velocity,
                                 Vector<MultiFab const*> const& density,
                                 Vector<MultiFab const*> const& tracer_old,
                                 Vector<MultiFab const*> const& tracer_new,
                                 bool include_pressure_gradient)
{
    for (int lev = 0; lev <= finest_level; ++lev)
       compute_vel_forces_on_level (lev, *vel_forces[lev], *velocity[lev], *density[lev],
                                         *tracer_old[lev], *tracer_new[lev], include_pressure_gradient);
}

void incflo::compute_vel_forces_on_level (int lev,
                                                MultiFab& vel_forces,
                                          const MultiFab& /*velocity*/,
                                          const MultiFab& density,
                                          const MultiFab& tracer_old,
                                          const MultiFab& tracer_new,
                                          bool include_pressure_gradient)
{
    GpuArray<Real,3> l_gravity{m_gravity[0],m_gravity[1],m_gravity[2]};
    GpuArray<Real,3> l_gp0{m_gp0[0], m_gp0[1], m_gp0[2]};

    auto const dx = geom[lev].CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
            Box const& bx = mfi.tilebox();
            Array4<Real>       const& vel_f =  vel_forces.array(mfi);
            Array4<Real const> const&   rho =     density.const_array(mfi);
            Array4<Real const> const& gradp = m_leveldata[lev]->gp.const_array(mfi);

            if (m_use_boussinesq) {
                // This uses a Boussinesq approximation where the buoyancy depends on
                //      first tracer rather than density
                Array4<Real const> const& tra_o = tracer_old.const_array(mfi);
                Array4<Real const> const& tra_n = tracer_new.const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    int n = 0; // Potential temperature

                    Real rhoinv = Real(1.0)/rho(i,j,k);
                    Real ft = Real(0.5) * (tra_o(i,j,k,n) + tra_n(i,j,k,n));

                    if (include_pressure_gradient)
                    {
                        AMREX_D_TERM(vel_f(i,j,k,0) = -gradp(i,j,k,0)*rhoinv + l_gravity[0] * ft;,
                                     vel_f(i,j,k,1) = -gradp(i,j,k,1)*rhoinv + l_gravity[1] * ft;,
                                     vel_f(i,j,k,2) = -gradp(i,j,k,2)*rhoinv + l_gravity[2] * ft;);
                    } else {
                        AMREX_D_TERM(vel_f(i,j,k,0) =                          l_gravity[0] * ft;,
                                     vel_f(i,j,k,1) =                          l_gravity[1] * ft;,
                                     vel_f(i,j,k,2) =                          l_gravity[2] * ft;);
                    }
                });

            } else if (m_probtype == 16) {
                Real Re = 1./m_mu;  // Note this assumes you are running exactly the problem set up, with U = 1 and L = 1 and rho = 1.
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhoinv = Real(1.0)/rho(i,j,k);

                    if (include_pressure_gradient)
                    {
                        AMREX_D_TERM(vel_f(i,j,k,0) = -(gradp(i,j,k,0)+l_gp0[0])*rhoinv + l_gravity[0];,
                                     vel_f(i,j,k,1) = -(gradp(i,j,k,1)+l_gp0[1])*rhoinv + l_gravity[1];,
                                     vel_f(i,j,k,2) = -(gradp(i,j,k,2)+l_gp0[2])*rhoinv + l_gravity[2];);
                    } else {
                        AMREX_D_TERM(vel_f(i,j,k,0) = -(               l_gp0[0])*rhoinv + l_gravity[0];,
                                     vel_f(i,j,k,1) = -(               l_gp0[1])*rhoinv + l_gravity[1];,
                                     vel_f(i,j,k,2) = -(               l_gp0[2])*rhoinv + l_gravity[2];);
                    }

                    Real x = (i+0.5) * dx[0];
                    Real y = (j+0.5) * dx[1];

                    Real f     = x*x*x*x - 2.*x*x*x + x*x;
                    Real g     = y*y*y*y - y*y;
                    Real capF  =  0.2 * x*x*x*x*x   -  0.5 * x*x*x*x   + (1./3.)* x*x*x;
                    Real capF1 = -4.0 * x*x*x*x*x*x + 12.0 * x*x*x*x*x - 14.    * x*x*x*x + 8.0 * x*x*x - 2.0 * x*x;
                    Real capF2 = 0.5 * f * f;
                    Real capG1 = -24.0 * y*y*y*y*y + 8.0 * y*y*y - 4.0 * y;

                    Real  fp   = 4.0 * x*x*x - 6.0*x*x + 2.0*x;
                    //Real  fpp  = 12.0 * x*x - 12.0*x + 2.0;
                    Real  fppp = 24.0 * x - 12.0;

                    Real  gp   =  4.0 * y*y*y - 2.0*y;
                    Real  gpp  = 12.0 * y*y - 2.0;

                    vel_f(i,j,k,1) += 8.0 / Re * (24.0 * capF + 2.0 * fp * gpp + fppp * g) + 64.0 * (capF2 * capG1 - g * gp * capF1);
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhoinv = Real(1.0)/rho(i,j,k);

                    if (include_pressure_gradient)
                    {
                        AMREX_D_TERM(vel_f(i,j,k,0) = -(gradp(i,j,k,0)+l_gp0[0])*rhoinv + l_gravity[0];,
                                     vel_f(i,j,k,1) = -(gradp(i,j,k,1)+l_gp0[1])*rhoinv + l_gravity[1];,
                                     vel_f(i,j,k,2) = -(gradp(i,j,k,2)+l_gp0[2])*rhoinv + l_gravity[2];);
                    } else {
                        AMREX_D_TERM(vel_f(i,j,k,0) = -(               l_gp0[0])*rhoinv + l_gravity[0];,
                                     vel_f(i,j,k,1) = -(               l_gp0[1])*rhoinv + l_gravity[1];,
                                     vel_f(i,j,k,2) = -(               l_gp0[2])*rhoinv + l_gravity[2];);
                    }
                });
            }
    }
///////////////////////////////////////////////////////////////////////////	
	// add surface tension
	// surface tension = sigma*kappa*grad(VOF)/rho
	// sigma: surface tension coefficient
	// kappa: curvature of the interface
	// grad(VOF): gradient of VOF field variable
	// rho:   density
	
	//fixme: we just consider the surface tension for first tracer	
    if (m_vof_advect_tracer && m_sigma[0]!=0.){
	
	  VolumeOfFluid*  vof_p = get_volume_of_fluid ();
	   
      //choice 1: The original cell-centered kappa and rho are averaged to face center. Grad(VOF) and  
	  // surface tension (SF) are calculated at face center. Then the face-centered SF is finally averaged to cell center. 
	  //choice 2: Similar to choice 1, SF is estimated at the face center and then averaged to the cell nodes.
	  // the node-centered SF is finally averaged to cell center.
	  //Choice 3: The original cell-centered rho and VOF are averaged to nodes. The node-centered rho is averaged
      // to cell center. grad(VOF) is estimated at the center of the cell edge and then averaged to the cell center.
      // finally, SF is calculated using original cell-centered kappa, averaged rho, and averaged grad(VOF).
      //Choice 4: SF is calculated using original cell-centered kappa, rho and center-difference for grad(VOF). 	  
	  
	  int choice = 3;
	  
	  
	  const auto& ba = density.boxArray();
      const auto& dm = density.DistributionMap();
      const auto& fact = density.Factory();
      Array<MultiFab,AMREX_SPACEDIM> face_val{AMREX_D_DECL(
	                                 MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                              dm, 1, 0, MFInfo(), fact))};
	  // store the nodal values (the last component stores the node-centered VOF)
      MultiFab node_val(amrex::convert(ba,IntVect::TheNodeVector()),dm, AMREX_SPACEDIM+1, 0 , MFInfo(), fact);		
      if (choice==1) {
    	// kappa, rho, and grad(VOF) are first averaged to the center of cell faces.
        // The cell-centered force is then obtained by averaging the face-centered value. 						  
	    average_cellcenter_to_face(GetArrOfPtrs(face_val), density, Geom(lev));
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          face_val[idim].invert(m_sigma[0], 0);
        }
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(vel_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {      
         // Note nodaltilebox will not include the nodal index beyond boundaries between neighboring
	     // titles. Therefore,if we want to use face values (i.e., face_val) immediately below (commented
	     // out), we must create index space for the face-centered values of the tiled region 
         // (i.e., surroundingNodes()). 	   
	     Box const& bx  = mfi.tilebox();
		 AMREX_D_TERM(
	       Box const& xbx = surroundingNodes(bx,0);,
	       Box const& ybx = surroundingNodes(bx,1);,
	       Box const& zbx = surroundingNodes(bx,2););	   
	     Array4<Real const> const& tra   = tracer_new.const_array(mfi);
         Array4<Real const> const& kap   = vof_p->kappa[lev].const_array(mfi);
	     AMREX_D_TERM(Array4<Real > const& xfv = face_val[0].array(mfi);,
			          Array4<Real > const& yfv = face_val[1].array(mfi);,
			          Array4<Real > const& zfv = face_val[2].array(mfi););
	     AMREX_D_TERM(
            ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              Real kaf;
		      if(kap(i,j,k,0)!=VOF_NODATA && kap(i-1,j,k,0)!=VOF_NODATA)
		  	    kaf=Real(0.5)*(kap(i,j,k,0)+kap(i-1,j,k,0));
		      else if (kap(i,j,k,0)!=VOF_NODATA)
		       kaf=kap(i,j,k);
		      else if (kap(i-1,j,k,0)!=VOF_NODATA)
               kaf=kap(i-1,j,k);
		      else
               kaf=0.;
		      xfv(i,j,k) *= kaf*(tra(i,j,k,0)-tra(i-1,j,k,0))/dx[0];				  		   			   
            });,      
	        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              Real kaf;
		      if(kap(i,j,k,0)!=VOF_NODATA && kap(i,j-1,k,0)!=VOF_NODATA)
		     	kaf=Real(0.5)*(kap(i,j,k,0)+kap(i,j-1,k,0));
		      else if (kap(i,j,k,0)!=VOF_NODATA)
		         kaf=kap(i,j,k);
		      else if (kap(i,j-1,k,0)!=VOF_NODATA)
                 kaf=kap(i,j-1,k);
		      else
                 kaf=0.;
		       yfv(i,j,k) *= kaf*(tra(i,j,k,0)-tra(i,j-1,k,0))/dx[1];				  		   			   
            });,	

	        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              Real kaf;
		      if(kap(i,j,k,0)!=VOF_NODATA && kap(i,j,k-1,0)!=VOF_NODATA)
		     	kaf=Real(0.5)*(kap(i,j,k,0)+kap(i,j,k-1,0));
		      else if (kap(i,j,k,0)!=VOF_NODATA)
		         kaf=kap(i,j,k);
		      else if (kap(i,j,k-1,0)!=VOF_NODATA)
                 kaf=kap(i,j,k-1);
		      else
                 kaf=0.;
		       zfv(i,j,k) *= kaf*(tra(i,j,k,0)-tra(i,j,k-1,0))/dx[2];	
	        
            });
          )	// end AMREX_D_TERM		
   
	     // immediately calculate the cell-centered surface tension force using the face-centered value 
	     // If we uncomment the following, we must use surroundingNodes() to define the index space (i.e.,
	     // xbx, ybx, zbx) for cell faces in the beginning. 
	      Array4<Real>  const& forarr  = vof_p->force[lev].array(mfi);
	      Array4<Real>  const& vel_f   = vel_forces.array(mfi);	   
          ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {	
		    		 
		    AMREX_D_TERM(vel_f(i,j,k,0) -= Real(0.5)*(xfv(i,j,k,0)+xfv(i+1,j,k,0));,
		   	          vel_f(i,j,k,1) -= Real(0.5)*(yfv(i,j,k,0)+yfv(i,j+1,k,0));,
		   	          vel_f(i,j,k,2) -= Real(0.5)*(zfv(i,j,k,0)+zfv(i,j,k+1,0)););	   
	        AMREX_D_TERM(forarr(i,j,k,0) = -Real(0.5)*(xfv(i,j,k,0)+xfv(i+1,j,k,0));,
		   	          forarr(i,j,k,1) = -Real(0.5)*(yfv(i,j,k,0)+yfv(i,j+1,k,0));,
		   	          forarr(i,j,k,2) = -Real(0.5)*(zfv(i,j,k,0)+zfv(i,j,k+1,0)););			
	      });
	    }
      } 
      else if (choice ==2){
       // kappa, rho, and grad(VOF) are first averaged to the center of cell faces.
       // The face-centered values are then averaged to the cell node.
       // Finally, the nodal values are averaged to the cell center.	   
	   average_cellcenter_to_face(GetArrOfPtrs(face_val), density, Geom(lev));
       for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
         face_val[idim].invert(m_sigma[0], 0);
       }
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(vel_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {      
       // Note nodaltilebox will not include the nodal index beyond boundaries between neighboring
	   // titles. Therefore,if we want to use face values (i.e., face_val) immediately below (commented
	   // out), we must create index space for the face-centered values of the tiled region 
       // (i.e., surroundingNodes()).Since we currently calculate all face values in all boxes and then 
       // convert them to node-centered value, it is better that we use (nodaltilebox()) to avoid the 
	   // repeated calculation of face-centered values at cell faces which are shared by two tiles.	   
	   Box const& xbx = mfi.nodaltilebox(0);
       Box const& ybx = mfi.nodaltilebox(1);
       Box const& zbx = mfi.nodaltilebox(2);	   
	   Array4<Real const> const& tra   = tracer_new.const_array(mfi);
       Array4<Real const> const& kap   = vof_p->kappa[lev].const_array(mfi);
	   AMREX_D_TERM( Array4<Real > const& xfv = face_val[0].array(mfi);,
			         Array4<Real > const& yfv = face_val[1].array(mfi);,
			         Array4<Real > const& zfv = face_val[2].array(mfi););
	   
	     AMREX_D_TERM(
            ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              Real kaf;
		      if(kap(i,j,k,0)!=VOF_NODATA && kap(i-1,j,k,0)!=VOF_NODATA)
		  	 kaf=Real(0.5)*(kap(i,j,k,0)+kap(i-1,j,k,0));
		      else if (kap(i,j,k,0)!=VOF_NODATA)
		       kaf=kap(i,j,k);
		      else if (kap(i-1,j,k,0)!=VOF_NODATA)
               kaf=kap(i-1,j,k);
		      else
               kaf=0.;
		      xfv(i,j,k) *= kaf*(tra(i,j,k,0)-tra(i-1,j,k,0))/dx[0];				  		   			   
            });,      
	        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              Real kaf;
		      if(kap(i,j,k,0)!=VOF_NODATA && kap(i,j-1,k,0)!=VOF_NODATA)
		     	kaf=Real(0.5)*(kap(i,j,k,0)+kap(i,j-1,k,0));
		      else if (kap(i,j,k,0)!=VOF_NODATA)
		         kaf=kap(i,j,k);
		      else if (kap(i,j-1,k,0)!=VOF_NODATA)
                 kaf=kap(i,j-1,k);
		      else
                 kaf=0.;
		       yfv(i,j,k) *= kaf*(tra(i,j,k,0)-tra(i,j-1,k,0))/dx[1];				  		   			   
            });,	
	        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              Real kaf;
		      if(kap(i,j,k,0)!=VOF_NODATA && kap(i,j,k-1,0)!=VOF_NODATA)
		     	kaf=Real(0.5)*(kap(i,j,k,0)+kap(i,j,k-1,0));
		      else if (kap(i,j,k,0)!=VOF_NODATA)
		         kaf=kap(i,j,k);
		      else if (kap(i,j,k-1,0)!=VOF_NODATA)
                 kaf=kap(i,j,k-1);
		      else
                 kaf=0.;
		       zfv(i,j,k) *= kaf*(tra(i,j,k,0)-tra(i,j,k-1,0))/dx[2];	
	        
            });
          )	// end AMREX_D_TERM	   
   
	   }
       static int oct[3][2] = { { 1, 2 }, { 0, 2 }, { 0, 1 } };
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(vel_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {        
	     Box const& bx  = mfi.tilebox();
	     //We will immediately use the nodal values so we must use surroundingNodes()instead of nodaltilebox()
	     //See the previous comments for calculating face-centered value.
	     Box const& nbx = surroundingNodes(bx);
	     Box const& vbx  = mfi.validbox();
		 AMREX_D_TERM(
	       Box const& xvbx = surroundingNodes(vbx,0);,
	       Box const& yvbx = surroundingNodes(vbx,1);,
	       Box const& zvbx = surroundingNodes(vbx,2););
	     AMREX_D_TERM( 
		   Array4<Real > const& xfv = face_val[0].array(mfi);,
		   Array4<Real > const& yfv = face_val[1].array(mfi);,
		   Array4<Real > const& zfv = face_val[2].array(mfi););
	     Array4<Real > const& nv = node_val.array(mfi);	   
         ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {			 
		    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim){
		      nv(i,j,k,dim)=0.;
		      int nt=0;
		      for (int detj = 0; detj < 2; ++detj) 
                for (int deti = 0; deti < 2; ++deti){
			      Array<int,3> in{i, j, k};  
			      in[oct[dim][0]]-=deti,in[oct[dim][1]]-=detj;
 		  	      if (dim==0&& xvbx.contains(in[0],in[1],in[2])){
		  	  	   nv(i,j,k,dim)+= xfv(in[0],in[1],in[2]); 
		  	  	   nt++;
		  	      }
		  	      else if (dim==1&& yvbx.contains(in[0],in[1],in[2])){
                       nv(i,j,k,dim)+= yfv(in[0],in[1],in[2]);	
		  	  	   nt++;
		  	      }				
		  	      else if (dim==2&& zvbx.contains(in[0],in[1],in[2])){
		  	  	   nv(i,j,k,dim)+= zfv(in[0],in[1],in[2]);
		  	  	   nt++;
		  	      }						   		      
			  }
		      if(nt>0) nv(i,j,k,dim)/= Real(nt);		   
	        }		 
	     });
	     Array4<Real>  const& forarr  = vof_p->force[lev].array(mfi);
	     Array4<Real>  const& vel_f   = vel_forces.array(mfi); 	   
         ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {									   
	  	    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim){			 
              vel_f(i,j,k,dim) -= Real(0.125)*(nv(i,j  ,k  ,dim) + nv(i+1,j  ,k  ,dim)
                                             + nv(i,j+1,k  ,dim) + nv(i+1,j+1,k  ,dim)
                                             + nv(i,j  ,k+1,dim) + nv(i+1,j  ,k+1,dim)
                                             + nv(i,j+1,k+1,dim) + nv(i+1,j+1,k+1,dim));
              forarr(i,j,k,dim) =-Real(0.125)*(nv(i,j  ,k  ,dim) + nv(i+1,j  ,k  ,dim)
                                             + nv(i,j+1,k  ,dim) + nv(i+1,j+1,k  ,dim)
                                             + nv(i,j  ,k+1,dim) + nv(i+1,j  ,k+1,dim)
                                             + nv(i,j+1,k+1,dim) + nv(i+1,j+1,k+1,dim));
	  	 								  
	  	    }			 		
	      });		 
       }	   	   
      } 
      else if (choice ==3){     	  
	   MultiFab center_val(ba,dm,2,0,MFInfo(), fact);
static int oct[3][2] = { { 1, 2 }, { 0, 2 }, { 0, 1 } };
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(vel_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {        
	     Box const& nbx = surroundingNodes(mfi.tilebox());
		 Box const& vbx = mfi.validbox();
	     Array4<Real > const& nv = node_val.array(mfi);
         Array4<Real const> const& tra   = tracer_new.const_array(mfi);	
         Array4<Real const> const& kappa = vof_p->kappa[lev].const_array(mfi);	
         Array4<Real const> const& rho   = density.const_array(mfi);	   
         ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {	
	       /*  if(i==8&&j==8&&k==8){ 
	         Print()<<"zbx   "<<"low  "<<tra(i,j,k,0)<<" high " <<tra(i,j,k-1,0)<<"  "
	         << "  "<<zfv(i,j,k)<<"\n";	     		 			 
		    }	*/	  		 
		 // calculate the node-centered VOF 
		   nv(i,j,k,AMREX_SPACEDIM)=0.;
		   nv(i,j,k,0)=0.;
		   nv(i,j,k,1)=0.;
		   int nt=0, nrho=0, nkap=0;		 
		   for (int detk = 0; detk > -2; --detk) 
            for (int detj = 0; detj > -2; --detj) 
             for (int deti = 0; deti > -2; --deti) {
		       Array<int,3> in{i+deti,j+detj,k+detk};
		       //if (vbx.contains(in[0],in[1],in[2])){
			   // averaging VOF to nodes
		  	   nv(i,j,k,AMREX_SPACEDIM)+= tra(in[0],in[1],in[2]); 
		  	   nt++;
			   //averaging kappa to nodes
		  	   if(kappa(in[0],in[1],in[2],0)!=VOF_NODATA){
		  		 nv(i,j,k,0)+= kappa(in[0],in[1],in[2],0); 
		  	     nkap++;  
		  	   }
			   //averaging density to nodes
		  	   nv(i,j,k,1)+= rho(in[0],in[1],in[2],0);
		  	   nrho++; 
		       //}			 
		     }
		   if(nt>0) nv(i,j,k,AMREX_SPACEDIM)/= Real(nt);
		   if(nkap>0)
			 nv(i,j,k,0)/= Real(nkap);
		   else 
			 nv(i,j,k,0)=0.;
		   nv(i,j,k,1)/= Real(nrho);		 
	     });
       }	   	   
	   average_node_to_cellcenter(center_val, 0, node_val, 0, 2);
       //MultiFab::Copy(density, center_val , 1, 0, 1, density.nGrow());	   
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(vel_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {        
	   Box const& bx  = mfi.tilebox();
	   //We will immediately use the nodal values so we must use surroundingNodes()instead of nodaltilebox()
	   //See the previous comments for calculating face-centered value.
	   Box const& nbx = surroundingNodes(bx);
	   Box const& vbx  = mfi.validbox();
	   Array4<Real > const& nv = node_val.array(mfi);
       Array4<Real const> const& tra   = tracer_new.const_array(mfi);	
       Array4<Real const> const& kappa = vof_p->kappa[lev].const_array(mfi);	
       Array4<Real const> const& rho   = density.const_array(mfi);		   
	   Array4<Real>  const& forarr  = vof_p->force[lev].array(mfi);
	   Array4<Real>  const& vel_f   = vel_forces.array(mfi); 
       Array4<Real>  const& center  = center_val.array(mfi);	   
       ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {	
		 Array<Real,AMREX_SPACEDIM> gradVof{0.};
		 for (int dim = 0; dim < AMREX_SPACEDIM; ++dim){          
		  for (int detj = 0; detj < 2; ++detj) 
           for (int deti = 0; deti < 2; ++deti){
			 Array<int,3> in0{i, j, k},in1{i, j, k};  
			 in0[oct[dim][0]]+=deti,in0[oct[dim][1]]+=detj;
             in1[oct[dim][0]]+=deti,in1[oct[dim][1]]+=detj;
             in1[dim] +=1;			 
		     gradVof[dim] +=Real(.25)*(nv(in1[0],in1[1],in1[2],AMREX_SPACEDIM)-
			                           nv(in0[0],in0[1],in0[2],AMREX_SPACEDIM))/dx[dim];						   
		   }
		   /*Array<int,3> in0{i, j, k},in1{i, j, k};
		   in1[dim] +=1,in0[dim] -=1;;
		   gradVof[dim] = Real(0.5)*(tra(in1[0],in1[1],in1[2],0)-tra(in0[0],in0[1],in0[2],0))/dx[dim];*/
         }		   
		 for (int dim = 0; dim < AMREX_SPACEDIM; ++dim){			 
           if(kappa(i,j,k,0)!=VOF_NODATA){
		     vel_f(i,j,k,dim) -= m_sigma[0]*kappa(i,j,k,0)/center(i,j,k,1)*gradVof[dim];
		     forarr(i,j,k,dim) =-m_sigma[0]*kappa(i,j,k,0)/center(i,j,k,1)*gradVof[dim];
		   }
           else {
             for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
              forarr(i,j,k,dim) = 0.;			  
		   }														  
		 }			 		
	   });
	
	
	}	
   }   
   else if(choice==4) {
//cell-centered surface tension force	
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(vel_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
          Box const& bx = mfi.tilebox();
          Array4<Real>       const& vel_f = vel_forces.array(mfi);
          Array4<Real const> const& rho   = density.const_array(mfi);
	      Array4<Real const> const& tra   = tracer_new.const_array(mfi);
          Array4<Real const> const& kappa = vof_p->kappa[lev].const_array(mfi);	
          Array4<Real >      const& forarr  = vof_p->force[lev].array(mfi);			  
          ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {		   
			if(kappa(i,j,k,0)!=VOF_NODATA){
              Real sig_kappa = m_sigma[0]*kappa(i,j,k,0)/rho(i,j,k);
             //note: the minus sign is used because of the way curvature is calculated			  
			  AMREX_D_TERM(
			   vel_f(i,j,k,0) -= Real(0.5)*(tra(i+1,j,k,0)-tra(i-1,j,k,0))/dx[0]*sig_kappa;,
			   vel_f(i,j,k,1) -= Real(0.5)*(tra(i,j+1,k,0)-tra(i,j-1,k,0))/dx[1]*sig_kappa;,
			   vel_f(i,j,k,2) -= Real(0.5)*(tra(i,j,k+1,0)-tra(i,j,k-1,0))/dx[2]*sig_kappa;);
			   
			  AMREX_D_TERM(
			   forarr(i,j,k,0) = -Real(0.5)*(tra(i+1,j,k,0)-tra(i-1,j,k,0))/dx[0]*sig_kappa;,
			   forarr(i,j,k,1) = -Real(0.5)*(tra(i,j+1,k,0)-tra(i,j-1,k,0))/dx[1]*sig_kappa;,
			   forarr(i,j,k,2) = -Real(0.5)*(tra(i,j,k+1,0)-tra(i,j,k-1,0))/dx[2]*sig_kappa;);			   
			   
			}
            else {
              for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
                forarr(i,j,k,dim) = 0.;				  
			}			
          });	   			   		
	    }

  }	
	
	
	
	
	}// end if (m_vof_advect_tracer)
}
