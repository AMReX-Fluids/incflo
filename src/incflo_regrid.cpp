#include <incflo.H>
#include <memory>

using namespace amrex;

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void incflo::MakeNewLevelFromCoarse (int lev,
                                     Real time,
                                     const BoxArray& ba,
                                     const DistributionMapping& dm)
{
    BL_PROFILE("incflo::MakeNewLevelFromCoarse()");

    if (m_verbose > 0) {
        amrex::Print() << "Making new level " << lev << " from coarse" << std::endl;
    }

#ifdef AMREX_USE_EB
    std::unique_ptr<FabFactory<FArrayBox> > new_fact = makeEBFabFactory(geom[lev], ba, dm,
                                                                        {nghost_eb_basic(),
                                                                         nghost_eb_volume(),
                                                                         nghost_eb_full()},
                                                                        EBSupport::full);
#else
    std::unique_ptr<FabFactory<FArrayBox> > new_fact(new FArrayBoxFactory());
#endif
    std::unique_ptr<LevelData> new_leveldata
        (new LevelData(ba, dm, *new_fact, this));
    fillcoarsepatch_velocity(lev, time, new_leveldata->velocity, 0);
    fillcoarsepatch_density(lev, time, new_leveldata->density, 0);
    if (m_ntrac > 0) {
        fillcoarsepatch_tracer(lev, time, new_leveldata->tracer, 0);
    }
    fillcoarsepatch_gradp(lev, time, new_leveldata->gp, 0);

    if (m_use_cc_proj) {
        new_leveldata->p_cc.setVal(0.0);
    } else {
        new_leveldata->p_nd.setVal(0.0);
    }

    m_leveldata[lev] = std::move(new_leveldata);
    m_factory[lev] = std::move(new_fact);

    m_diffusion_tensor_op.reset();
    m_diffusion_scalar_op.reset();
    if (m_vof_advect_tracer){
      std::unique_ptr<VolumeOfFluid::LevelData> new_leveldata_vof
                 (new VolumeOfFluid::LevelData(ba, dm, *m_factory[lev], this));
      get_volume_of_fluid()->m_leveldata[lev] = std::move(new_leveldata_vof);
    }
    //ptr_VOF.reset();

    // Note: finest_level has not yet been updated and so we use lev
#ifdef AMREX_USE_EB
    macproj = std::make_unique<Hydro::MacProjector>(Geom(0,lev),
                      MLMG::Location::FaceCentroid,  // Location of mac_vec
                      MLMG::Location::FaceCentroid,  // Location of beta
                      MLMG::Location::CellCenter  ); // Location of solution variable phi
#else
    macproj = std::make_unique<Hydro::MacProjector>(Geom(0,lev));
#endif
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void incflo::RemakeLevel (int lev, Real time, const BoxArray& ba,
                          const DistributionMapping& dm)
{
    BL_PROFILE("incflo::RemakeLevel()");

    if (m_verbose > 0) {
        amrex::Print() << "Remaking level " << lev << std::endl;
    }

#ifdef AMREX_USE_EB
    std::unique_ptr<FabFactory<FArrayBox> > new_fact = makeEBFabFactory(geom[lev], ba, dm,
                                                                        {nghost_eb_basic(),
                                                                         nghost_eb_volume(),
                                                                         nghost_eb_full()},
                                                                        EBSupport::full);
#else
    std::unique_ptr<FabFactory<FArrayBox> > new_fact(new FArrayBoxFactory());
#endif
    std::unique_ptr<LevelData> new_leveldata
        (new LevelData(ba, dm, *new_fact, this));
    fillpatch_velocity(lev, time, new_leveldata->velocity, 0);
    fillpatch_density(lev, time, new_leveldata->density, 0);
    if (m_ntrac > 0) {
// Updated to include the patch for one layer of ghost cells in the multifab 'tracer'.
// This is necessary because 'tracer_vof_advection' is called before 'tracer'
// is updated in incflo::ApplyPredictor. The VOF advection process requires
// ghost cell values to compute fluxes for cells adjacent to the boundary.
        fillpatch_tracer(lev, time, new_leveldata->tracer, 1);
    }

    fillpatch_gradp(lev, time, new_leveldata->gp, 0);

    if (m_use_cc_proj) {
        new_leveldata->p_cc.setVal(0.0);
    } else {
        new_leveldata->p_nd.setVal(0.0);
    }

    m_leveldata[lev] = std::move(new_leveldata);
    m_factory[lev] = std::move(new_fact);


//    for (MFIter mfi(m_leveldata[lev]->tracer); mfi.isValid(); ++mfi) {
//       Box const& bx = mfi./*growntilebox(1);*/validbox();
//       Array4<Real const> const& vof_arr = m_leveldata[lev]->tracer.const_array(mfi);
//       ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//       {
//         auto fvol = vof_arr(i,j,k,0);
//    Print()<<"regrid--lev= "<<lev<<"(i,j) "<<i<<","<<j<<bx;
//    Print()<<" vof= "<<fvol<<" \n";
//
//       }); //end ParallelFor
//    } //end MFIter

    //make_mixedBC_mask(lev, ba, dm);

    m_diffusion_tensor_op.reset();
    m_diffusion_scalar_op.reset();
    if (m_vof_advect_tracer){
      std::unique_ptr<VolumeOfFluid::LevelData> new_leveldata_vof
                 (new VolumeOfFluid::LevelData(ba, dm, *m_factory[lev], this));
      get_volume_of_fluid()->m_leveldata[lev] = std::move(new_leveldata_vof);
      //fixme: it may be better to move the following to incflo::Evolve() after calling regrid()
      auto& ldvof=*ptr_VOF->m_leveldata[lev]; /*VOF data for level lev*/
      ptr_VOF->tracer_vof_update(lev, m_leveldata[lev]->tracer, ldvof.height);
      ptr_VOF->curvature_calculation(lev, m_leveldata[lev]->tracer, ldvof.height, ldvof.kappa);
      //diffuse the VOF by averaging
      const auto& ba = m_leveldata[lev]->tracer.boxArray();
      const auto& dm = m_leveldata[lev]->tracer.DistributionMap();
      MultiFab tracer_df(ba,dm,1,m_leveldata[lev]->tracer.nGrow(),MFInfo(), *m_factory[lev]);
      MultiFab::Copy(tracer_df, m_leveldata[lev]->tracer, 0, 0, 1, m_leveldata[lev]->tracer.nGrow());
      for (int i=0;i<m_number_of_averaging;i++){
         ptr_VOF->variable_filtered(lev, tracer_df);
      }
      update_vof_density (lev, m_leveldata[lev]->density, tracer_df);

      //auto tag_vector_ptrs = ptr_VOF->get_vector_ptr([](VolumeOfFluid::LevelData& ld) -> MultiFab& {return ld.tag;});
      //ptr_VOF->domain_tag_droplets (lev, grids, geom, get_tracer_new(),tag_vector_ptrs);
    }
    //ptr_VOF.reset();

#ifdef AMREX_USE_EB
    macproj = std::make_unique<Hydro::MacProjector>(Geom(0,finest_level),
                      MLMG::Location::FaceCentroid,  // Location of mac_vec
                      MLMG::Location::FaceCentroid,  // Location of beta
                      MLMG::Location::CellCenter  ); // Location of solution variable phi
#else
    macproj = std::make_unique<Hydro::MacProjector>(Geom(0,finest_level));
#endif

#ifdef INCFLO_USE_PARTICLES
    particleData.Redistribute();
#endif
}

// Delete level data
// overrides the pure virtual function in AmrCore
void incflo::ClearLevel (int lev)
{
    BL_PROFILE("incflo::ClearLevel()");
    m_leveldata[lev].reset();
    m_factory[lev].reset();
    m_diffusion_tensor_op.reset();
    m_diffusion_scalar_op.reset();
    macproj.reset();
    if (m_vof_advect_tracer){
        ptr_VOF->m_leveldata[lev].reset();
    }
    //ptr_VOF.reset();

}
