#include <incflo.H>
#include <prob_bc.H>
#include <AMReX_FillPatchUtil.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBInterpolater.H>
#endif

using namespace amrex;

void incflo::fillpatch_velocity (int lev, Real time, MultiFab& vel, int ng)
{
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > physbc
            (geom[lev], get_velocity_bcrec(),
             IncfloVelFill{m_probtype, m_bc_velocity});
        FillPatchSingleLevel(vel, IntVect(ng), time,
                             {&(m_leveldata[lev]->velocity_o),
                              &(m_leveldata[lev]->velocity)},
                             {m_t_old[lev], m_t_new[lev]}, 0, 0, AMREX_SPACEDIM, geom[lev],
                             physbc, 0);
    } else {
        const auto& bcrec = get_velocity_bcrec();
        PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > cphysbc
            (geom[lev-1], bcrec, IncfloVelFill{m_probtype, m_bc_velocity});
        PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > fphysbc
            (geom[lev], bcrec, IncfloVelFill{m_probtype, m_bc_velocity});
#ifdef AMREX_USE_EB
        Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
        Interpolater* mapper = &cell_cons_interp;
#endif
        switch (m_fillpatch_method){
          case 0:
            // This FillPatch operation does not use the ghost cells of the coarser level
            FillPatchTwoLevels(vel, IntVect(ng), time,
                               {&(m_leveldata[lev-1]->velocity_o),
                                &(m_leveldata[lev-1]->velocity)},
                               {m_t_old[lev-1], m_t_new[lev-1]},
                               {&(m_leveldata[lev]->velocity_o),
                                &(m_leveldata[lev]->velocity)},
                               {m_t_old[lev], m_t_new[lev]},
                               0, 0, AMREX_SPACEDIM, geom[lev-1], geom[lev],
                               cphysbc, 0, fphysbc, 0,
                               refRatio(lev-1), mapper, bcrec, 0);
            break;
          case 1:
            // This FillPatch operation interpolates using the ghost cells of the coarser level
            // via `PhysBCFunctUseCoarseGhost`, which is defined in `AMReX_PhysBCFunct.h`.
            // For implementation details, see `AMReX_FillPatchUtil_I.h`.
            //
            // When the `blocking_factor` is small (e.g., 1, 2, or 4), specifically used for generating
            // quad-/octree-like grids, this FillPatch method is necessary instead of the previous one.
            FillPatchTwoLevels (vel, IntVect(ng), IntVect (0), time,
                                {&(m_leveldata[lev-1]->velocity_o),
                                 &(m_leveldata[lev-1]->velocity)},
                                {m_t_old[lev-1], m_t_new[lev-1]},
                                {&(m_leveldata[lev]->velocity_o),
                                 &(m_leveldata[lev]->velocity)},
                                {m_t_old[lev], m_t_new[lev]},
                                0, 0, AMREX_SPACEDIM, geom[lev-1], geom[lev],
                                refRatio(lev-1), mapper, bcrec, 0);
            //The physical boundary condition is not enforced in the above fillpatch, so we have to do it here.
            fphysbc.FillBoundary(vel, 0, AMREX_SPACEDIM, IntVect(ng), time, 0);
            break;
          case 2:
            //for quad-/Oct-tree like grids, it is safter to use FillPatchNLevels
            Vector<PhysBCFunct<GpuBndryFuncFab<IncfloVelFill>>> physbcs;
            for (int ilev = 0; ilev <= finest_level; ++ilev) {
                physbcs.emplace_back(geom[ilev],bcrec,IncfloVelFill{m_probtype, m_bc_velocity});
            }
            Vector<Vector<MultiFab*>> smf(finest_level+1);
            Vector<Vector<Real>> st(finest_level+1);
            for (int ilev = 0; ilev <= finest_level; ++ilev) {
              smf[ilev] = {&(m_leveldata[ilev]->velocity_o), &(m_leveldata[ilev]->velocity)};
              st[ilev] = {m_t_old[ilev], m_t_new[ilev]};
            }
            FillPatchNLevels(vel, lev, IntVect(ng), time, smf, st, 0, 0, AMREX_SPACEDIM, geom,
                             physbcs, 0, ref_ratio, mapper, bcrec, 0);
            break;
        }
    }
}

void incflo::fillpatch_density (int lev, Real time, MultiFab& density, int ng)
{
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > physbc(geom[lev], get_density_bcrec(),
                                                            IncfloDenFill{m_probtype, m_bc_density, m_bc_velocity});
        FillPatchSingleLevel(density, IntVect(ng), time,
                             {&(m_leveldata[lev]->density_o),
                              &(m_leveldata[lev]->density)},
                             {m_t_old[lev], m_t_new[lev]}, 0, 0, 1, geom[lev],
                             physbc, 0);
    } else {
        const auto& bcrec = get_density_bcrec();
        PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > cphysbc
            (geom[lev-1], bcrec, IncfloDenFill{m_probtype, m_bc_density, m_bc_velocity});
        PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > fphysbc
            (geom[lev], bcrec, IncfloDenFill{m_probtype, m_bc_density, m_bc_velocity});
#ifdef AMREX_USE_EB
        Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
        Interpolater* mapper = &cell_cons_interp;
#endif

        switch (m_fillpatch_method){
          case 0:
            // This FillPatch operation does not use the ghost cells of the coarser level
            FillPatchTwoLevels(density, IntVect(ng), time,
                               {&(m_leveldata[lev-1]->density_o),
                                &(m_leveldata[lev-1]->density)},
                               {m_t_old[lev-1], m_t_new[lev-1]},
                               {&(m_leveldata[lev]->density_o),
                                &(m_leveldata[lev]->density)},
                               {m_t_old[lev], m_t_new[lev]},
                               0, 0, 1, geom[lev-1], geom[lev],
                               cphysbc, 0, fphysbc, 0,
                               refRatio(lev-1), mapper, bcrec, 0);
            break;
          case 1:
            // This FillPatch operation interpolates using the ghost cells of the coarser level
            // via `PhysBCFunctUseCoarseGhost`, which is defined in `AMReX_PhysBCFunct.h`.
            // For implementation details, see `AMReX_FillPatchUtil_I.h`.
            //
            // When the `blocking_factor` is small (e.g., 1, 2, or 4), specifically used for generating
            // quad-/octree-like grids, this FillPatch method is necessary instead of the previous one.
            FillPatchTwoLevels (density, IntVect(ng), IntVect (0), time,
                                {&(m_leveldata[lev-1]->density_o),
                                 &(m_leveldata[lev-1]->density)},
                                {m_t_old[lev-1], m_t_new[lev-1]},
                                {&(m_leveldata[lev]->density_o),
                                 &(m_leveldata[lev]->density)},
                                {m_t_old[lev], m_t_new[lev]},
                                0, 0, 1, geom[lev-1], geom[lev],
                                refRatio(lev-1), mapper, bcrec, 0);
            //The physical boundary condition is not enforced in the above fillpatch, so we have to do it here.
            fphysbc.FillBoundary(density, 0, 1, IntVect(ng), time, 0);
            break;
          case 2:
            //for quad-/Oct-tree like grids, it is safter to use FillPatchNLevels
            Vector<PhysBCFunct<GpuBndryFuncFab<IncfloDenFill>>> physbcs;
            for (int ilev = 0; ilev <= finest_level; ++ilev) {
                physbcs.emplace_back(geom[ilev],bcrec,IncfloDenFill{m_probtype, m_bc_density, m_bc_velocity});
            }
            Vector<Vector<MultiFab*>> smf(finest_level+1);
            Vector<Vector<Real>> st(finest_level+1);
            for (int ilev = 0; ilev <= finest_level; ++ilev) {
              smf[ilev] = {&(m_leveldata[ilev]->density_o), &(m_leveldata[ilev]->density)};
              st[ilev] = {m_t_old[ilev], m_t_new[ilev]};
            }
            FillPatchNLevels(density, lev, IntVect(ng), time, smf, st, 0, 0, 1, geom,
                             physbcs, 0, ref_ratio, mapper, bcrec, 0);
            break;
        }
    }
}

void incflo::fillpatch_tracer (int lev, Real time, MultiFab& tracer, int ng)
{
    if (m_ntrac <= 0) return;
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloTracFill> > physbc
            (geom[lev], get_tracer_bcrec(), IncfloTracFill{m_probtype, m_ntrac, m_bc_tracer_d, m_bc_velocity});
        FillPatchSingleLevel(tracer, IntVect(ng), time,
                             {&(m_leveldata[lev]->tracer_o),
                              &(m_leveldata[lev]->tracer)},
                             {m_t_old[lev], m_t_new[lev]}, 0, 0, m_ntrac, geom[lev],
                             physbc, 0);
    } else {
        const auto& bcrec = get_tracer_bcrec();
        PhysBCFunct<GpuBndryFuncFab<IncfloTracFill> > cphysbc
            (geom[lev-1], bcrec, IncfloTracFill{m_probtype, m_ntrac, m_bc_tracer_d, m_bc_velocity});
        PhysBCFunct<GpuBndryFuncFab<IncfloTracFill> > fphysbc
            (geom[lev], bcrec, IncfloTracFill{m_probtype, m_ntrac, m_bc_tracer_d, m_bc_velocity});
#ifdef AMREX_USE_EB
        Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
        Interpolater* mapper = &cell_cons_interp;
#endif
        switch (m_fillpatch_method){
          case 0:
            // This FillPatch operation does not use the ghost cells of the coarser level
            FillPatchTwoLevels(tracer, IntVect(ng), time,
                             {&(m_leveldata[lev-1]->tracer_o),
                              &(m_leveldata[lev-1]->tracer)},
                             {m_t_old[lev-1], m_t_new[lev-1]},
                             {&(m_leveldata[lev]->tracer_o),
                              &(m_leveldata[lev]->tracer)},
                             {m_t_old[lev], m_t_new[lev]},
                             0, 0, m_ntrac, geom[lev-1], geom[lev],
                             cphysbc, 0, fphysbc, 0,
                             refRatio(lev-1), mapper, bcrec, 0);
            break;
          case 1:
            // This FillPatch operation interpolates using the ghost cells of the coarser level
            // via `PhysBCFunctUseCoarseGhost`, which is defined in `AMReX_PhysBCFunct.h`.
            // For implementation details, see `AMReX_FillPatchUtil_I.h`.
            //
            // When the `blocking_factor` is small (e.g., 1, 2, or 4), specifically used for generating
            // quad-/octree-like grids, this FillPatch method is necessary instead of the previous one.
            FillPatchTwoLevels (tracer, IntVect(ng), IntVect (0), time,
                                {&(m_leveldata[lev-1]->tracer_o),
                                 &(m_leveldata[lev-1]->tracer)},
                                {m_t_old[lev-1], m_t_new[lev-1]},
                                {&(m_leveldata[lev]->tracer_o),
                                 &(m_leveldata[lev]->tracer)},
                                {m_t_old[lev], m_t_new[lev]},
                                0, 0, m_ntrac, geom[lev-1], geom[lev],
                                refRatio(lev-1), mapper, bcrec, 0);
            //The physical boundary condition is not enforced in the above fillpatch, so we have to do it here.
            fphysbc.FillBoundary(tracer, 0, m_ntrac, IntVect(ng), time, 0);
            break;
          case 2:
            //for quad-/Oct-tree like grids, it is safter to use FillPatchNLevels
            // FillPatchNLevels is not as fast as fillpatch operation in case 1
            if (m_vof_advect_tracer){
              mapper = &(get_volume_of_fluid()->vof_interp);
              get_volume_of_fluid()->vof_interp.lev = lev;
            }
            Vector<BCRec> tmp_bcrec;
            if (m_vof_advect_tracer){
              tmp_bcrec.reserve(bcrec.size() + AMREX_SPACEDIM+1);
              tmp_bcrec.insert(tmp_bcrec.end(), bcrec.begin(), bcrec.end());
              const auto& bcrec_force = get_force_bcrec();
              tmp_bcrec.insert(tmp_bcrec.end(), bcrec_force.begin(), bcrec_force.begin() + AMREX_SPACEDIM+1);
            }
            else{
                tmp_bcrec=bcrec;
            }
            Vector<PhysBCFunct<GpuBndryFuncFab<IncfloTracFill>>> physbcs;
            for (int ilev = 0; ilev <= finest_level; ++ilev) {
              physbcs.emplace_back(geom[ilev],tmp_bcrec,IncfloTracFill{m_probtype, m_ntrac, m_bc_tracer_d, m_bc_velocity});
            }
            Vector<Vector<MultiFab*>> smf(finest_level+1);
            Vector<Vector<Real>> st(finest_level+1);
            if (m_vof_advect_tracer){
               Vector<MultiFab> vof_all(finest_level+1);

               for (int ilev = 0; ilev <= finest_level; ++ilev) {
                 vof_all[ilev].define(grids[ilev], dmap[ilev], 2+AMREX_SPACEDIM, nghost_state(), MFInfo());
                 Copy(vof_all[ilev], m_leveldata[ilev]->tracer, 0, 0, 1, ng);
                 Copy(vof_all[ilev], ptr_VOF->m_leveldata[ilev]->normal, 0, 1, AMREX_SPACEDIM, ng);
                 Copy(vof_all[ilev], ptr_VOF->m_leveldata[ilev]->alpha, 0, AMREX_SPACEDIM+1, 1, ng);
                 smf[ilev] = {&(vof_all[ilev])};
                 st[ilev] = {time};
               //smf[ilev] = {&(m_leveldata[ilev]->tracer_o), &(m_leveldata[ilev]->tracer)};
              //st[ilev] = {m_t_old[ilev], m_t_new[ilev]};
               }
               //note  tracer_temp is a temporary multifab created using the 'tracer's BoxArray and DistributionMappings
               //which may be different than 'm_leveldata[lev]->tracer'. For example, a new grid is generated in RemakeLevel()
               //before  fillpatch() is called.
               MultiFab tracer_tmp(tracer.boxArray(), tracer.DistributionMap(), 2+AMREX_SPACEDIM, tracer.nGrow(), MFInfo());
               FillPatchNLevels(tracer_tmp, lev, IntVect(ng), time, smf, st, 0, 0, AMREX_SPACEDIM+2, geom,
                                physbcs, 0, ref_ratio, mapper, tmp_bcrec, 0);
               Copy(tracer, tracer_tmp, 0, 0, 1, ng);
               //we only need to copy the data of current level lev.
               // Note: The 'tracer' and 'vof_all[lev]' may be created using different BoxArrays and DistributionMappings.
               // For example, when fillpatch_tracer() is called after a new grid is generated in RemakeLevel(),
               // 'tracer' may have a different BoxArray than 'm_leveldata[lev]->tracer'.
               // Therefore, ParallelCopy() is used here instead of Copy() to ensure proper data transfer.
               //tracer.ParallelCopy(tracer_tmp, 0, 0, 1, nghost_state(), nghost_state(), geom[lev].periodicity());
//    if (lev == 2 && m_nstep==31)
//    for (MFIter mfi(tracer); mfi.isValid(); ++mfi) {
//       Box const& bx = mfi./*growntilebox(ng);*/validbox();
//       Array4<Real const> const& vof_arr = tracer.const_array(mfi);
//       //Array4<Real const> const& vof_arr0 = tracer.const_array(mfi);
//       ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//       {
//         auto fvol = vof_arr(i,j,k,0);
//    Print()<<"fp-lev= "<<lev<<"(i,j) "<<i<<","<<j<<bx;
//    Print()<<" vof= "<<fvol/*<<"vof_copy= "<<vof_arr0(i,j,k,0)*/<<" \n";
//
//       }); //end ParallelFor
//    } //end MFIter



           }
            else{
              for (int ilev = 0; ilev <= finest_level; ++ilev) {
                 smf[ilev] = {&(m_leveldata[ilev]->tracer)};
                 st[ilev] = {time};
              }
              FillPatchNLevels(tracer, lev, IntVect(ng), time, smf, st, 0, 0, m_ntrac, geom,
                               physbcs, 0, ref_ratio, mapper, bcrec, 0);
            }
            break;
        }
    }
}

void incflo::fillpatch_gradp (int lev, Real time, MultiFab& gp, int ng)
{
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloForFill> > physbc
            (geom[lev], get_force_bcrec(), IncfloForFill{m_probtype});
        FillPatchSingleLevel(gp, IntVect(ng), time,
                             {&(m_leveldata[lev]->gp)}, {time},
                             0, 0, AMREX_SPACEDIM, geom[lev], physbc, 0);
    } else {
        const auto& bcrec = get_force_bcrec();
        PhysBCFunct<GpuBndryFuncFab<IncfloForFill> > cphysbc
            (geom[lev-1], bcrec, IncfloForFill{m_probtype});
        PhysBCFunct<GpuBndryFuncFab<IncfloForFill> > fphysbc
            (geom[lev], bcrec, IncfloForFill{m_probtype});
#ifdef AMREX_USE_EB
        Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
        Interpolater* mapper = &cell_cons_interp;
#endif

        switch (m_fillpatch_method){
          case 0:
            // This FillPatch operation does not use the ghost cells of the coarser level
            FillPatchTwoLevels(gp, IntVect(ng), time,
                               {&(m_leveldata[lev-1]->gp)}, {time},
                               {&(m_leveldata[lev]->gp)}, {time},
                               0, 0, AMREX_SPACEDIM, geom[lev-1], geom[lev],
                               cphysbc, 0, fphysbc, 0,
                               refRatio(lev-1), mapper, bcrec, 0);
            break;
          case 1:
            // This FillPatch operation interpolates using the ghost cells of the coarser level
            // via `PhysBCFunctUseCoarseGhost`, which is defined in `AMReX_PhysBCFunct.h`.
            // For implementation details, see `AMReX_FillPatchUtil_I.h`.
            //
            // When the `blocking_factor` is small (e.g., 1, 2, or 4), specifically used for generating
            // quad-/octree-like grids, this FillPatch method is necessary instead of the previous one.
            FillPatchTwoLevels (gp, IntVect(ng), IntVect (0), time,
                                {&(m_leveldata[lev-1]->gp)}, {time},
                                {&(m_leveldata[lev]->gp)}, {time},
                                0, 0, AMREX_SPACEDIM, geom[lev-1], geom[lev],
                                refRatio(lev-1), mapper, bcrec, 0);
            //The physical boundary condition is not enforced in the above fillpatch, so we have to do it here.
            fphysbc.FillBoundary(gp, 0, AMREX_SPACEDIM, IntVect(ng), time, 0);
            break;
          case 2:
            //for quad-/Oct-tree like grids, it is safter to use FillPatchNLevels
            Vector<PhysBCFunct<GpuBndryFuncFab<IncfloForFill>>> physbcs;
            for (int ilev = 0; ilev <= finest_level; ++ilev) {
                physbcs.emplace_back(geom[ilev],bcrec,IncfloForFill{m_probtype});
            }
            Vector<Vector<MultiFab*>> smf(finest_level+1);
            Vector<Vector<Real>> st(finest_level+1);
            for (int ilev = 0; ilev <= finest_level; ++ilev) {
              smf[ilev] = {&(m_leveldata[ilev]->gp)};
              st[ilev] = {time};
            }
            FillPatchNLevels(gp, lev, IntVect(ng), time, smf, st, 0, 0, AMREX_SPACEDIM, geom,
                             physbcs, 0, ref_ratio, mapper, bcrec, 0);
            break;
        }

    }
}

void incflo::fillpatch_force (Real time, Vector<MultiFab*> const& force, int ng)
{
    const int ncomp = force[0]->nComp();
    const auto& bcrec = get_force_bcrec();
    int lev = 0;
    {
        PhysBCFunct<GpuBndryFuncFab<IncfloForFill> > physbc
            (geom[lev], bcrec, IncfloForFill{m_probtype});
        FillPatchSingleLevel(*force[lev], IntVect(ng), time,
                             {force[lev]}, {time},
                             0, 0, ncomp, geom[lev],
                             physbc, 0);
    }
    for (lev = 1; lev <= finest_level; ++lev)
    {
        PhysBCFunct<GpuBndryFuncFab<IncfloForFill> > cphysbc
            (geom[lev-1], bcrec, IncfloForFill{m_probtype});
        PhysBCFunct<GpuBndryFuncFab<IncfloForFill> > fphysbc
            (geom[lev  ], bcrec, IncfloForFill{m_probtype});
        Interpolater* mapper = &pc_interp;

        switch (m_fillpatch_method){
          case 0:
            // This FillPatch operation does not use the ghost cells of the coarser level
            FillPatchTwoLevels(*force[lev], IntVect(ng), time,
                               {force[lev-1]}, {time},
                               {force[lev  ]}, {time},
                               0, 0, ncomp, geom[lev-1], geom[lev],
                               cphysbc, 0, fphysbc, 0,
                               refRatio(lev-1), mapper, bcrec, 0);
            break;
          case 1:
            // This FillPatch operation interpolates using the ghost cells of the coarser level
            // via `PhysBCFunctUseCoarseGhost`, which is defined in `AMReX_PhysBCFunct.h`.
            // For implementation details, see `AMReX_FillPatchUtil_I.h`.
            //
            // When the `blocking_factor` is small (e.g., 1, 2, or 4), specifically used for generating
            // quad-/octree-like grids, this FillPatch method is necessary instead of the previous one.
            FillPatchTwoLevels(*force[lev], IntVect(ng), IntVect (0), time,
                               {force[lev-1]}, {time},
                               {force[lev  ]}, {time},
                                0, 0, ncomp, geom[lev-1], geom[lev],
                                refRatio(lev-1), mapper, bcrec, 0);
            //The physical boundary condition is not enforced in the above fillpatch, so we have to do it here.
            fphysbc.FillBoundary(*force[lev], 0, ncomp, IntVect(ng), time, 0);
            break;
          case 2:
            //for quad-/Oct-tree like grids, it is safter to use FillPatchNLevels
            Vector<PhysBCFunct<GpuBndryFuncFab<IncfloForFill>>> physbcs;
            for (int ilev = 0; ilev <= finest_level; ++ilev) {
                physbcs.emplace_back(geom[ilev],bcrec,IncfloForFill{m_probtype});
            }
            Vector<Vector<MultiFab*>> smf(finest_level+1);
            Vector<Vector<Real>> st(finest_level+1);
            for (int ilev = 0; ilev <= finest_level; ++ilev) {
              smf[ilev] = {force[lev-1]};
              st[ilev] = {time};
            }
            FillPatchNLevels(*force[lev], lev, IntVect(ng), time, smf, st, 0, 0, ncomp, geom,
                             physbcs, 0, ref_ratio, mapper, bcrec, 0);
            break;
        }
    }
}

void incflo::fillcoarsepatch_velocity (int lev, Real time, MultiFab& vel, int ng)
{
    const auto& bcrec = get_velocity_bcrec();
    PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > cphysbc
        (geom[lev-1], bcrec, IncfloVelFill{m_probtype, m_bc_velocity});
    PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > fphysbc
        (geom[lev], bcrec, IncfloVelFill{m_probtype, m_bc_velocity});
#ifdef AMREX_USE_EB
    Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
        (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
    Interpolater* mapper = &cell_cons_interp;
#endif
    amrex::InterpFromCoarseLevel(vel, IntVect(ng), time,
                                 m_leveldata[lev-1]->velocity, 0, 0, AMREX_SPACEDIM,
                                 geom[lev-1], geom[lev],
                                 cphysbc, 0, fphysbc, 0,
                                 refRatio(lev-1), mapper, bcrec, 0);
}

void incflo::fillcoarsepatch_density (int lev, Real time, MultiFab& density, int ng)
{
    const auto& bcrec = get_density_bcrec();
    PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > cphysbc
        (geom[lev-1], bcrec, IncfloDenFill{m_probtype, m_bc_density, m_bc_velocity});
    PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > fphysbc
        (geom[lev], bcrec, IncfloDenFill{m_probtype, m_bc_density, m_bc_velocity});
#ifdef AMREX_USE_EB
    Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
        (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
    Interpolater* mapper = &cell_cons_interp;
#endif
    amrex::InterpFromCoarseLevel(density, IntVect(ng), time,
                                 m_leveldata[lev-1]->density, 0, 0, 1,
                                 geom[lev-1], geom[lev],
                                 cphysbc, 0, fphysbc, 0,
                                 refRatio(lev-1), mapper, bcrec, 0);
}

void incflo::fillcoarsepatch_tracer (int lev, Real time, MultiFab& tracer, int ng)
{
    if (m_ntrac <= 0) return;

    const auto& bcrec = get_tracer_bcrec();
    PhysBCFunct<GpuBndryFuncFab<IncfloTracFill> > cphysbc
        (geom[lev-1], bcrec, IncfloTracFill{m_probtype, m_ntrac, m_bc_tracer_d, m_bc_velocity});
    PhysBCFunct<GpuBndryFuncFab<IncfloTracFill> > fphysbc
        (geom[lev], bcrec, IncfloTracFill{m_probtype, m_ntrac, m_bc_tracer_d, m_bc_velocity});
#ifdef AMREX_USE_EB
    Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
        (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
    Interpolater* mapper = &cell_cons_interp;
#endif
    amrex::InterpFromCoarseLevel(tracer, IntVect(ng), time,
                                 m_leveldata[lev-1]->tracer, 0, 0, m_ntrac,
                                 geom[lev-1], geom[lev],
                                 cphysbc, 0, fphysbc, 0,
                                 refRatio(lev-1), mapper, bcrec, 0);
}

void incflo::fillcoarsepatch_gradp (int lev, Real time, MultiFab& gp, int ng)
{
    const auto& bcrec = get_force_bcrec();
    PhysBCFunct<GpuBndryFuncFab<IncfloForFill> > cphysbc
        (geom[lev-1], bcrec, IncfloForFill{m_probtype});
    PhysBCFunct<GpuBndryFuncFab<IncfloForFill> > fphysbc
        (geom[lev], bcrec, IncfloForFill{m_probtype});
#ifdef AMREX_USE_EB
    Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
        (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
    Interpolater* mapper = &cell_cons_interp;
#endif
    amrex::InterpFromCoarseLevel(gp, IntVect(ng), time,
                                 m_leveldata[lev-1]->gp, 0, 0, AMREX_SPACEDIM,
                                 geom[lev-1], geom[lev],
                                 cphysbc, 0, fphysbc, 0,
                                 refRatio(lev-1), mapper, bcrec, 0);
}
