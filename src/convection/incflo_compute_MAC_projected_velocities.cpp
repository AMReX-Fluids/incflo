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

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        mac_phi[lev]->FillBoundary(geom[lev].periodicity());

#ifdef AMREX_USE_EB
        const EBFArrayBoxFactory* ebfact = &EBFactory(lev);
#endif

        // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
        //    returned from this call are on face CENTROIDS
        HydroUtils::ExtrapVelToFaces(*vel[lev], *vel_forces[lev],
                                      AMREX_D_DECL(*u_mac[lev], *v_mac[lev], *w_mac[lev]),
                                      get_velocity_bcrec(), get_velocity_bcrec_device_ptr(),
                                      geom[lev], l_dt,
#ifdef AMREX_USE_EB
                                      *ebfact,
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

    if (m_flow_through_eb) {
       for (int lev=0; lev <= finest_level; ++lev)
       {
          macproj->setEBDirichlet(lev, *get_velocity_eb()[lev]);
       }
    }

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
