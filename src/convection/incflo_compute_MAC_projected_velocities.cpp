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
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel_f(i,j,k,0) += divtau(i,j,k,0)/rho(i,j,k);,
                                         vel_f(i,j,k,1) += divtau(i,j,k,1)/rho(i,j,k);,
                                         vel_f(i,j,k,2) += divtau(i,j,k,2)/rho(i,j,k););
                        });
                    }
                    else {
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel_f(i,j,k,0) += divtau(i,j,k,0);,
                                         vel_f(i,j,k,1) += divtau(i,j,k,1);,
                                         vel_f(i,j,k,2) += divtau(i,j,k,2););
                        });
                    } // end m_advect_momentum

                } // end MFIter

            } // end lev

        } // end m_godunov_include_diff_in_forcing

        if (nghost_force() > 0)
            fillpatch_force(m_cur_time, vel_forces, nghost_force());

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
        amrex::average_cellcenter_to_face(GetArrOfPtrs(inv_rho[lev]), *density[lev], geom[lev]);
#endif

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            inv_rho[lev][idim].invert(l_dt, 0);
        }
    }

    //
    // Initialize (or redefine the beta in) the MacProjector
    //
    Print()<<"Mac proj needs init? "<<macproj->needInitialization()<<std::endl;
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

        // Not sure exactly where this goes. Does it only need initialization? set every time???
        // initProj only take as many levels as 1/rho, so if a level is added, would have to
        // re-init proj anyway.
        // robin_a, etc needs to exist to finest_level and that's it
        // Values need to be in the ghost cells, although the BC is considered on face.
        if ( m_has_mixedBC ) {
            int nghost = 1;

            for (int lev = 0; lev <= finest_level; ++lev)
            {
                // define and fill Robin BC a, b, and f
                MultiFab robin_a (grids[lev],dmap[lev],1,nghost,MFInfo(),Factory(lev));
                MultiFab robin_b (grids[lev],dmap[lev],1,nghost,MFInfo(),Factory(lev));
                MultiFab robin_f (grids[lev],dmap[lev],1,nghost,MFInfo(),Factory(lev));

                // fixme - for vis, init robin MFs, although i don't know that this is generally needed
                robin_a = 0;
                robin_b = 0;
                robin_f = 0;

                // I don't believe interior values get used,
                // only ghost cells of robin BC arrays are used in MLMG
                // bc in ghost cells that are outside the domain.
                // What does MLMG expect for periodic? Should have been taken
                // care of in creating the mask
                //amrex::Geometry const& gm = Geom(lev);
                //Box const& domain = gm.growPeriodicDomain(nghost);
                Box const& domain = Geom(lev).Domain();
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    Orientation olo(dir,Orientation::low);
                    Orientation ohi(dir,Orientation::high);

                    if (m_bc_type[olo] == BC::mixed || m_bc_type[ohi] == BC::mixed) {
                        Box dlo = (m_bc_type[olo] == BC::mixed) ? adjCellLo(domain,dir) : Box();
                        Box dhi = (m_bc_type[ohi] == BC::mixed) ? adjCellHi(domain,dir) : Box();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                    // FIXME - don't think we want tiling here...
                        for (MFIter mfi(robin_a,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                            Box const& gbx = amrex::grow(mfi.validbox(),nghost);
                            Box blo = gbx & dlo;
                            Box bhi = gbx & dhi;
                            Array4<Real> const& a_arr     = robin_a.array(mfi);
                            Array4<Real> const& b_arr     = robin_b.array(mfi);
                            Array4<Real> const& f_arr     = robin_f.array(mfi);
                            //Array4<int const> const& mask = m_mixedBC_mask[lev]->const_array(mfi);

                            // Robin BC:   a u + b du/dn = f  -- inflow,  Neumann   a=0, b=1, f=0
                            //                                -- outflow, Dirichlet a=1, b=0, f=0
// Only need to fill Robin BC sides, MLMG will check for Robin BC first
                            if (blo.ok()) {
                                // robin_b is the same as the nodalBC_mask, except with vals in
                                // cell-centered ghosts rather than on BC face.
                                // So, lo side i_cc->i_nd+1, hi side i_cc=i_nd
                                //Dim3 shift = IntVect::TheDimensionVector(dir).dim3();

                                prob_set_MAC_robinBCs(olo, blo, a_arr, b_arr, f_arr, lev);

                                // amrex::ParallelFor(blo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                                // {
                                //     //b_arr(i,j,k) = mask(i+shift.x, j+shift.y, k+shift.z);
                                //     //robin_a is the "opposite" of b; 0->1 and vice versa
                                //     a_arr(i,j,k) = (mask(i+shift.x, j+shift.y, k+shift.z) == 1) ? 0. : 1.;
                                //     f_arr(i,j,k) = 0.;
                                // });
                            }
                            if (bhi.ok()) {

                                prob_set_MAC_robinBCs(ohi, bhi, a_arr, b_arr, f_arr, lev);

                                // amrex::ParallelFor(bhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                                // {
                                //     b_arr(i,j,k) = mask(i,j,k);
                                //     a_arr(i,j,k) = (mask(i,j,k) == 1) ? 0. : 1.;
                                //     f_arr(i,j,k) = 0.;
                                // });
                            }
                        }
                    }
                }

                VisMF::Write(robin_a, "ra");
                VisMF::Write(robin_b, "rb");
                VisMF::Write(robin_f, "rf");
                macproj->setLevelBC(lev, nullptr, &robin_a, &robin_b, &robin_f);
            }
        }
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
        if ( l_advection_type == "BDS" ) {
            l_advection_type = "Godunov";
        }

// FIXME - need to create the BC MF here. This is a iMF vs a Real MF, so can't use
        // the robin BC info from MAC. This is because if we want this to be extensible to
        // full position dependent BCs, would need the MF to hold the BCType, which are enums,
// so ints...
        //
        std::unique_ptr<iMultiFab> BC_MF;
        if (m_has_mixedBC) {
            BC_MF = make_BC_MF(lev, m_bcrec_velocity_d, "velocity");
        }
//fixme
    VisMF::Write(*vel[0],"vin");
    amrex::Write(*BC_MF,"bcmf");
//

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
                                      l_advection_type, PPM::default_limiter, BC_MF.get());
    }

//fixme
    VisMF::Write(*u_mac[0],"umac");
    VisMF::Write(*v_mac[0],"vmac");
//
    Vector<Array<MultiFab*,AMREX_SPACEDIM> > mac_vec(finest_level+1);
    for (int lev=0; lev <= finest_level; ++lev)
    {
        AMREX_D_TERM(mac_vec[lev][0] = u_mac[lev];,
                     mac_vec[lev][1] = v_mac[lev];,
                     mac_vec[lev][2] = w_mac[lev];);
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
            //mac_phi[lev]->mult(0.5,0,1,1);

        macproj->project(mac_phi,m_mac_mg_rtol,m_mac_mg_atol);

        for (int lev=0; lev <= finest_level; ++lev)
            mac_phi[lev]->mult(2.0,0,1,1);
    } else {
        macproj->project(m_mac_mg_rtol,m_mac_mg_atol);
    }
//fixme
    VisMF::Write(*u_mac[0],"umacp");
    VisMF::Write(*v_mac[0],"vmacp");
    // Note that the macproj->project call above ensures that the MAC velocities are averaged down --
    //      we don't need to do that again here
}
