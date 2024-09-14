#include <hydro_utils.H>
#include <hydro_bcs_K.H>
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

        if (m_godunov_include_diff_in_forcing) {

            for (int lev = 0; lev <= finest_level; ++lev) {
                auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*density[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box const& bx = mfi.tilebox();
                    Array4<Real> const& vel_f          = vel_forces[lev]->array(mfi);
                    Array4<Real const> const& rho      = density[lev]->array(mfi);
                    Array4<Real const> const& divtau   = ld.divtau_o.const_array(mfi);
                    if (m_advect_momentum) {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel_f(i,j,k,0) += divtau(i,j,k,0)/rho(i,j,k);,
                                         vel_f(i,j,k,1) += divtau(i,j,k,1)/rho(i,j,k);,
                                         vel_f(i,j,k,2) += divtau(i,j,k,2)/rho(i,j,k););
                        });
                    }
                    else {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel_f(i,j,k,0) += divtau(i,j,k,0);,
                                         vel_f(i,j,k,1) += divtau(i,j,k,1);,
                                         vel_f(i,j,k,2) += divtau(i,j,k,2););
                        });
                    } // end m_advect_momentum

                } // end MFIter

            } // end lev

        } // end m_godunov_include_diff_in_forcing

        if (nghost_force() > 0) {
            fillpatch_force(m_cur_time, vel_forces, nghost_force());
        }

    } // end m_advection_type


    // This will hold (1/rho) on faces
    Vector<Array<MultiFab,AMREX_SPACEDIM> > inv_rho(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(
           inv_rho[lev][0].define(u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),Factory(lev));,
           inv_rho[lev][1].define(v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),Factory(lev));,
           inv_rho[lev][2].define(w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),Factory(lev)));
    }

    for (int lev=0; lev <= finest_level; ++lev)
    {
#ifdef AMREX_USE_EB
        EB_interp_CellCentroid_to_FaceCentroid (*density[lev], GetArrOfPtrs(inv_rho[lev]),
                                                0, 0, 1, geom[lev], get_density_bcrec());
#else
        average_cellcenter_to_face(GetArrOfPtrs(inv_rho[lev]), *density[lev], geom[lev]);
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
        macproj->setDomainBC(get_mac_projection_bc(Orientation::low), get_mac_projection_bc(Orientation::high));

        if ( m_has_mixedBC ) {
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                auto const robin = make_robinBC_MFs(lev);
                macproj->setLevelBC(lev, nullptr,
                                    &robin[0], &robin[1], &robin[2]);
            }
        }
    } else {
#ifndef AMREX_USE_EB
        if (m_constant_density&&!m_vof_advect_tracer) {
            macproj->updateBeta(l_dt/m_ro_0);  // unnecessary unless m_ro_0 changes.
        } else
#endif
        {
            macproj->updateCoeffs(GetVecOfArrOfConstPtrs(inv_rho));
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        mac_phi[lev]->FillBoundary(geom[lev].periodicity());

#ifdef AMREX_USE_EB
        const EBFArrayBoxFactory* ebfact = &EBFactory(lev);
#endif

        //
        // For now, BDS is only for edge state prediction given a known advective velocity,
        // not this velocity extrapolation. So we'll use Godunov here instead.
        //
        std::string l_advection_type = m_advection_type;
        if ( m_advection_type == "BDS" ) {
            l_advection_type = "Godunov";
        }

        std::unique_ptr<iMultiFab> BC_MF;
        if (m_has_mixedBC) {
            // Create a MF to hold the BCType info. Note that this is different than the
            // bcs for the MAC projection because the MAC operates on phi, this is velocity.
            //
            // TODO? Could consider creating an incflo member variable to save the BC_MF
            //     amrex::Vector<std::unique_ptr<amrex::iMultiFab> > m_BC_MF;
            // Could stack components as vel, density, tracer, and then could use for
            // scalars' advective step as well. But not sure it really matters one way or
            // the other.
            BC_MF = make_BC_MF(lev, m_bcrec_velocity_d, "velocity");
        }

        // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
        //    returned from this call are on face CENTROIDS
        bool allow_inflow_on_outflow = false;
        HydroUtils::ExtrapVelToFaces(*vel[lev], *vel_forces[lev],
                                      AMREX_D_DECL(*u_mac[lev], *v_mac[lev], *w_mac[lev]),
                                      get_velocity_bcrec(), get_velocity_bcrec_device_ptr(),
                                      geom[lev], l_dt,
#ifdef AMREX_USE_EB
                                      *ebfact,
                                      m_eb_flow.enabled ? get_velocity_eb()[lev] : nullptr,
#endif
                                      m_godunov_ppm, m_godunov_use_forces_in_trans,
                                      l_advection_type, PPM::default_limiter,
                                      allow_inflow_on_outflow, BC_MF.get());

        //add surface tension
        //if(m_vof_advect_tracer)
        //  get_volume_of_fluid ()->velocity_face_source(lev,l_dt, AMREX_D_DECL(*u_mac[lev], *v_mac[lev], *w_mac[lev]));

    }

    Vector<Array<MultiFab*,AMREX_SPACEDIM> > mac_vec(finest_level+1);
    for (int lev=0; lev <= finest_level; ++lev)
    {
        AMREX_D_TERM(mac_vec[lev][0] = u_mac[lev];,
                     mac_vec[lev][1] = v_mac[lev];,
                     mac_vec[lev][2] = w_mac[lev];);
    }

    const auto *bc_vel_d = get_velocity_bcrec_device_ptr();

    std::unique_ptr<iMultiFab> velBC_MF;

    //
    // This loop is only necessary if there is time-dependent inflow but we don't have a good test for that
    //      so we do it if there is any inflow face
    //
    for (int lev=0; lev <= finest_level; ++lev)
    {
        MultiFab time_dep_inflow_vel(vel[lev]->boxArray(),vel[lev]->DistributionMap(),AMREX_SPACEDIM,1);
        time_dep_inflow_vel.setVal(0.);
        fillphysbc_velocity(lev, m_cur_time+0.5*l_dt, time_dep_inflow_vel, 1);

        Box domain(geom[lev].Domain());
        const auto dlo = lbound(domain);
        const auto dhi = ubound(domain);

        if (m_has_mixedBC) {
            velBC_MF = make_BC_MF(lev, m_bcrec_velocity_d, "velocity");
        }

        for (MFIter mfi(time_dep_inflow_vel,false); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.validbox();
            auto   cc_arr = time_dep_inflow_vel.const_array(mfi);

            Array4<int const> const& velbc_arr = velBC_MF ? (*velBC_MF).const_array(mfi)
                                                          : Array4<int const>{};

            auto umac_arr = mac_vec[lev][0]->array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int n = 0;
                const auto bc = HydroBC::getBC(i, j, k, n, domain, bc_vel_d, velbc_arr);
                if (i == dlo.x && ( bc.lo(0) == BCType::ext_dir ||
                                   (bc.lo(0) == BCType::direction_dependent && cc_arr(i-1,j,k,0) >= Real(0.0)) ) ) {
                    umac_arr(i,j,k) = cc_arr(i-1,j,k,0);
                }
                if (i == dhi.x && ( bc.hi(0) == BCType::ext_dir ||
                                   (bc.hi(0) == BCType::direction_dependent && cc_arr(i+1,j,k,0) <= Real(0.0)) ) ) {
                    umac_arr(i+1,j,k) = cc_arr(i+1,j,k,0);
                }
            });
            auto vmac_arr = mac_vec[lev][1]->array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int n = 1;
                const auto bc = HydroBC::getBC(i, j, k, n, domain, bc_vel_d, velbc_arr);
                if (j == dlo.y && ( bc.lo(1) == BCType::ext_dir||
                                   (bc.lo(1) == BCType::direction_dependent && cc_arr(i,j-1,k,1) >= Real(0.0)) ) ) {
                    vmac_arr(i,j,k) = cc_arr(i,j-1,k,1);
                }
                if (j == dhi.y && ( bc.hi(1) == BCType::ext_dir ||
                                   (bc.hi(1) == BCType::direction_dependent && cc_arr(i,j+1,k,1) <= Real(0.0)) ) ) {
                    vmac_arr(i,j+1,k) = cc_arr(i,j+1,k,1);
                }
            });
#if (AMREX_SPACEDIM > 2)
            auto wmac_arr = mac_vec[lev][2]->array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int n = 2;
                const auto bc = HydroBC::getBC(i, j, k, n, domain, bc_vel_d, velbc_arr);
                if (k == dlo.z && ( bc.lo(2) == BCType::ext_dir ||
                                   (bc.lo(2) == BCType::direction_dependent && cc_arr(1,j,k-1,2) >= Real(0.0)) ) ) {
                    wmac_arr(i,j,k) = cc_arr(i,j,k-1,2);
                }
                if (k == dhi.z && ( bc.hi(2) == BCType::ext_dir ||
                                   (bc.hi(2) == BCType::direction_dependent && cc_arr(i,j,k+1,2) <= Real(0.0)) ) ) {
                    wmac_arr(i,j,k+1) = cc_arr(i,j,k+1,2);
                }
            });
#endif
        } // mfi
    } // lev

    // Enforce solvability by matching outflow to inflow.
    if (has_inout_bndry)
    {
        HydroUtils::enforceInOutSolvability(mac_vec, get_velocity_bcrec().data(), geom);
    }

    macproj->setUMAC(mac_vec);

#ifdef AMREX_USE_EB
    if (m_eb_flow.enabled) {
       for (int lev=0; lev <= finest_level; ++lev)
       {
          macproj->setEBInflowVelocity(lev, *get_velocity_eb()[lev]);
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

        macproj->project(mac_phi,m_mac_mg_rtol,m_mac_mg_atol);

        for (int lev=0; lev <= finest_level; ++lev)
            mac_phi[lev]->mult(2.0,0,1,1);
    } else {
        macproj->project(m_mac_mg_rtol,m_mac_mg_atol);
    }
    // Note that the macproj->project call above ensures that the MAC velocities are averaged down --
    //      we don't need to do that again here
}
