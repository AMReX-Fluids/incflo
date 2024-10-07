#include <incflo.H>

using namespace amrex;

void diffuse_density(int lev, Vector<MultiFab*> const& density)
{
    const auto& ba = density[lev]->boxArray();
    const auto& dm = density[lev]->DistributionMap();
    const auto& fact = density[lev]->Factory();
    // store the nodal values (the last component stores the node-centered VOF)
    MultiFab node_val(amrex::convert(ba,IntVect::TheNodeVector()),dm, 1, 0 , MFInfo(), fact);
    MultiFab center_val(ba,dm,1,0,MFInfo(), fact);
    return;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(*density[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
         Box const& nbx = surroundingNodes(mfi.tilebox());
         Box const& vbx = mfi.validbox();
         Array4<Real > const& nv = node_val.array(mfi);
         Array4<Real const> const& rho   = density[lev]->const_array(mfi);
         ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
           nv(i,j,k,0)=0.;
           int nt=0, nrho=0, detk;
#if AMREX_SPACEDIM==3
           for (detk = 0; detk > -2; --detk)
#endif
            for (int detj = 0; detj > -2; --detj)
             for (int deti = 0; deti > -2; --deti) {
               Array<int,3> in{i+deti,j+detj,k+detk};
               //averaging density to nodes
                 nv(i,j,k,0)+= rho(in[0],in[1],in[2],0);
                 nrho++;
             }
           nv(i,j,k,1)/= Real(nrho);
         });
       }
       //average_node_to_cellcenter(center_val, 0, node_val, 0, 0);
      // Copy(*density[lev], center_val , 0, 0, 1, 0);

}

void incflo::update_density (StepType step_type)
{
    BL_PROFILE("incflo::update_density");

    int ng = (step_type == StepType::Corrector) ? 0 : 1;

    Real l_dt = m_dt;

    if (!m_constant_density)
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
          if(m_update_density_from_vof){
           update_vof_density (lev, get_density_new(),get_tracer_new());
           diffuse_density(lev, get_density_new());
          }
          else{
            auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real const> const& rho_old  = ld.density_o.const_array(mfi);
                Array4<Real      > const& rho_new  = ld.density.array(mfi);
                Array4<Real const> const& drdt_o = ld.conv_density_o.const_array(mfi);

                if (step_type == StepType::Predictor) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        rho_new(i,j,k) = rho_old(i,j,k) + l_dt * drdt_o(i,j,k);
                    });

                } else if (step_type == StepType::Corrector) {
                    Array4<Real const> const& drdt_n = ld.conv_density.const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        rho_new(i,j,k) = rho_old(i,j,k) + l_dt * Real(0.5) * (drdt_o(i,j,k) + drdt_n(i,j,k));
                    });
                }
            } // mfi
          }
        } // lev

        // Average down solution
        for (int lev = finest_level-1; lev >= 0; --lev) {
#ifdef AMREX_USE_EB
            amrex::EB_average_down(m_leveldata[lev+1]->density, m_leveldata[lev]->density,
                                   0, 1, refRatio(lev));
#else
            amrex::average_down(m_leveldata[lev+1]->density, m_leveldata[lev]->density,
                                0, 1, refRatio(lev));
#endif
        }

        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];

            // Fill ghost cells of new-time density if needed (we assume ghost cells of old density are already filled)
            if (ng > 0) {
                fillpatch_density(lev, m_t_new[lev], ld.density, ng);
            }

            // Define half-time density after the average down
            MultiFab::LinComb(ld.density_nph, Real(0.5), ld.density, 0, Real(0.5), ld.density_o, 0, 0, 1, ng);
        }

    } else {
        for (int lev = 0; lev <= finest_level; lev++) {
          if (m_vof_advect_tracer){
              //when VOF method is used to advect the tracer, density and viscosity of each cell will
              //depend the VOF field value of the cell.
              //fixme
             update_vof_density (lev, get_density_new(),get_tracer_old());
             //MultiFab::Copy(m_leveldata[lev]->density, m_leveldata[lev]->density_o, 0, 0, 1, m_leveldata[lev]->density_o.nGrow());
          }
          MultiFab::Copy(m_leveldata[lev]->density_nph, m_leveldata[lev]->density_o, 0, 0, 1, ng);
        }
    }
}
