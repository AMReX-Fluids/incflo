#include <hydro_utils.H>
#include <incflo.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void
incflo::compute_MAC_projected_velocities (
                                 Vector<MultiFab const*> const& vel,
                                 Vector<MultiFab const*> const& density,
                                 AMREX_D_DECL(Vector<MultiFab*> const& u_mac,
                                              Vector<MultiFab*> const& v_mac,
                                              Vector<MultiFab*> const& w_mac),
                                 Vector<MultiFab*> const& vel_forces,
                                 Real /*time*/)
{
    BL_PROFILE("incflo::compute_MAC_projected_velocities()");
    Real l_dt = m_dt;

    auto mac_phi = get_mac_phi();

    // We first compute the velocity forcing terms to be used in predicting
    //    to faces before the MAC projection
    if (m_advection_type != "MOL") {

        if (m_godunov_include_diff_in_forcing)
            for (int lev = 0; lev <= finest_level; ++lev)
                MultiFab::Add(*vel_forces[lev], m_leveldata[lev]->divtau_o, 0, 0, AMREX_SPACEDIM, 0);

        if (nghost_force() > 0)
            fillpatch_force(m_cur_time, vel_forces, nghost_force());
    }


    // This will hold (1/rho) on faces
    Vector<Array<MultiFab,AMREX_SPACEDIM> > inv_rho(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(
           inv_rho[lev][0].define(u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),u_mac[lev]->Factory());,
           inv_rho[lev][1].define(v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),v_mac[lev]->Factory());,
           inv_rho[lev][2].define(w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),w_mac[lev]->Factory()));
    }

    for (int lev=0; lev <= finest_level; ++lev)
    {
#ifdef AMREX_USE_EB
        EB_interp_CellCentroid_to_FaceCentroid (*density[lev], GetArrOfPtrs(inv_rho[lev]),
                                                0, 0, 1, geom[lev], get_density_bcrec());
#else
        amrex::average_cellcenter_to_face(GetArrOfPtrs(inv_rho[lev]), *density[lev], geom[lev]);
#endif

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            inv_rho[lev][idim].invert(l_dt, 0);
        }
    }

    //
    // Initialize (or redefine the beta in) the MacProjector
    //
    if (macproj->needInitialization())
    {
        LPInfo lp_info;
        lp_info.setMaxCoarseningLevel(m_mac_mg_max_coarsening_level);
#ifndef AMREX_USE_EB
        if (m_constant_density) {
            Vector<BoxArray> ba;
            Vector<DistributionMapping> dm;
            for (auto const& ir : inv_rho) {
                ba.push_back(ir[0].boxArray());
                dm.push_back(ir[0].DistributionMap());
            }
            macproj->initProjector(ba, dm, lp_info, l_dt/m_ro_0);
        } else
#endif
        {
            macproj->initProjector(lp_info, GetVecOfArrOfConstPtrs(inv_rho));
        }
        macproj->setDomainBC(get_projection_bc(Orientation::low), get_projection_bc(Orientation::high));
    } else {
#ifndef AMREX_USE_EB
        if (m_constant_density) {
            macproj->updateBeta(l_dt/m_ro_0);  // unnecessary unless m_ro_0 changes.
        } else
#endif
        {
            macproj->updateBeta(GetVecOfArrOfConstPtrs(inv_rho));
        }
    }

    // try to print out inv_rho
    auto const& irhox = inv_rho[0][0].const_array(0);
    auto const& irhoy = inv_rho[0][1].const_array(0);
   
    for (int i = 4; i < 13; i++){
        for (int j = 4; j < 13; j++){
            amrex::Print() << "irho" << IntVect(i,j) << ": " << irhox(i,j,0) << ", " << irhoy(i,j,0) << std::endl;
        }
        }

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        mac_phi[lev]->FillBoundary(geom[lev].periodicity());

#ifdef AMREX_USE_EB
        const EBFArrayBoxFactory* ebfact = &OldEBFactory(lev);
#endif

        // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
        //    returned from this call are on face CENTROIDS
        HydroUtils::ExtrapVelToFaces(*vel[lev], *vel_forces[lev],
                                      AMREX_D_DECL(*u_mac[lev], *v_mac[lev], *w_mac[lev]),
                                      get_velocity_bcrec(), get_velocity_bcrec_device_ptr(),
                                      geom[lev], l_dt,
#ifdef AMREX_USE_EB
                                      *ebfact,
                                      m_eb_flow.enabled ? get_velocity_eb()[lev] : nullptr,
#endif
                                      m_godunov_ppm, m_godunov_use_forces_in_trans,
                                      m_advection_type);
    }

    Vector<Array<MultiFab*,AMREX_SPACEDIM> > mac_vec(finest_level+1);
    for (int lev=0; lev <= finest_level; ++lev)
    {
        AMREX_D_TERM(mac_vec[lev][0] = u_mac[lev];,
                     mac_vec[lev][1] = v_mac[lev];,
                     mac_vec[lev][2] = w_mac[lev];);
    }

    macproj->setUMAC(mac_vec);

#ifdef AMREX_USE_EB
    Vector<MultiFab*> eb_vel_mod;

    if (m_eb_flow.enabled) {
       for (int lev=0; lev <= finest_level; ++lev)
       {
          MultiFab* new_mf = new MultiFab(grids[lev],dmap[lev],AMREX_SPACEDIM,0);
          eb_vel_mod.push_back(new_mf);
          MultiFab::Copy(*eb_vel_mod[lev], *get_velocity_eb()[lev], 0, 0, AMREX_SPACEDIM, 0);

          auto const& vfrac_old  = OldEBFactory(lev).getVolFrac();
          auto const& vfrac_new  =    EBFactory(lev).getVolFrac();
          auto const& bndry_area = OldEBFactory(lev).getBndryArea();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(*eb_vel_mod[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
              Box const& bx = mfi.tilebox();

              // Face-centered areas
              AMREX_D_TERM(const auto& apx = OldEBFactory(lev).getAreaFrac()[0]->const_array(mfi);,
                           const auto& apy = OldEBFactory(lev).getAreaFrac()[1]->const_array(mfi);,
                           const auto& apz = OldEBFactory(lev).getAreaFrac()[2]->const_array(mfi););
 
              Array4<Real      > const&   vel_arr  = eb_vel_mod[lev]->array(mfi);
              Array4<Real const> const& vfold_arr  =  vfrac_old.const_array(mfi);
              Array4<Real const> const& vfnew_arr  =  vfrac_new.const_array(mfi);
              Array4<Real const> const& ebarea_arr =  bndry_area.const_array(mfi);

              Real dx = geom[lev].CellSize()[0]; 

              amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                  if (vfold_arr(i,j,k) > 0. && vfold_arr(i,j,k) < 1.0) { 

                      AMREX_D_TERM(const Real dapx = apx(i+1,j  ,k  ) - apx(i,j,k);,
                                   const Real dapy = apy(i  ,j+1,k  ) - apy(i,j,k);,
                                   const Real dapz = apz(i  ,j  ,k+1) - apz(i,j,k););
#if (AMREX_SPACEDIM == 2)
                      Real apnorm = std::sqrt(dapx*dapx+dapy*dapy);
#elif (AMREX_SPACEDIM == 3)
                      Real apnorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
#endif
                      Real apnorm_inv = 1.0/apnorm;

                      AMREX_D_TERM(Real nx = dapx * apnorm_inv;,
                                   Real ny = dapy * apnorm_inv;,
                                   Real nz = dapz * apnorm_inv;);

                      Real delta_vol = vfnew_arr(i,j,k) - vfold_arr(i,j,k);
                     
                      AMREX_D_TERM(vel_arr(i,j,k,0) = -dx * delta_vol / l_dt / ebarea_arr(i,j,k) * nx;,
                                   vel_arr(i,j,k,1) = -dx * delta_vol / l_dt / ebarea_arr(i,j,k) * ny;,
                                   vel_arr(i,j,k,2) = -dx * delta_vol / l_dt / ebarea_arr(i,j,k) * nz;);
                  }
              });
          }

          //macproj->setEBInflowVelocity(lev, *get_velocity_eb()[lev]);
          macproj->setEBInflowVelocity(lev, *eb_vel_mod[lev]);
       }
    }
#endif

    if (m_verbose > 0) amrex::Print() << "MAC Projection:\n";
    //
    // Perform MAC projection
    //
    if (m_use_mac_phi_in_godunov)
    {
        // The MAC projection always starts with phi == 0, but we might like
        //     to change that so we reduce the cost of the MAC projection
        for (int lev=0; lev <= finest_level; ++lev)
            mac_phi[lev]->setVal(0.);
            //mac_phi[lev]->mult(0.5,0,1,1);

        macproj->project(mac_phi,m_mac_mg_rtol,m_mac_mg_atol);

        for (int lev=0; lev <= finest_level; ++lev)
            mac_phi[lev]->mult(2.0,0,1,1);
    } else {
        macproj->project(m_mac_mg_rtol,m_mac_mg_atol);
    }

    // Note that the macproj->project call above ensures that the MAC velocities are averaged down --
    //      we don't need to do that again here
}
