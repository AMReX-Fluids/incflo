#include <Godunov.H>
#ifdef AMREX_USE_EB
#include <EBGodunov.H>
#include <Redistribution.H>
#endif
#include <MOL.H>
#include <incflo.H>

using namespace amrex;

void incflo::init_advection ()
{
#ifdef AMREX_USE_EB
    if (m_advection_type == "Godunov")
    {
       // We default to conservative if doing Godunov
       m_iconserv_velocity.resize(AMREX_SPACEDIM, 1);
       m_iconserv_velocity_d.resize(AMREX_SPACEDIM, 1);
    } else {
       // We default to non-conservative if doing MOL
       m_iconserv_velocity.resize(AMREX_SPACEDIM, 0);
       m_iconserv_velocity_d.resize(AMREX_SPACEDIM, 0);
    }
#else
   m_iconserv_velocity.resize(AMREX_SPACEDIM, 0);
   m_iconserv_velocity_d.resize(AMREX_SPACEDIM, 0);
#endif

    m_iconserv_density.resize(1, 1);
    m_iconserv_density_d.resize(1, 1);

    m_iconserv_tracer.resize(m_ntrac, 1);
    m_iconserv_tracer_d.resize(m_ntrac, 1);
}

void
incflo::compute_convective_term (Vector<MultiFab*> const& conv_u,
                                 Vector<MultiFab*> const& conv_r,
                                 Vector<MultiFab*> const& conv_t,
                                 Vector<MultiFab const*> const& vel,
                                 Vector<MultiFab const*> const& density,
                                 Vector<MultiFab const*> const& tracer,
                                 AMREX_D_DECL(Vector<MultiFab*> const& u_mac,
                                              Vector<MultiFab*> const& v_mac,
                                              Vector<MultiFab*> const& w_mac),
                                 Vector<MultiFab*      > const& vel_forces,
                                 Vector<MultiFab*      > const& tra_forces,
                                 Real time)
{
    int ngmac = nghost_mac();

    // We first compute the velocity forcing terms to be used in predicting
    //    to faces before the MAC projection
    if (m_advection_type != "MOL") {

        bool include_pressure_gradient = !(m_use_mac_phi_in_godunov);
        compute_vel_forces(vel_forces, vel, density, tracer, tracer, include_pressure_gradient);

        if (m_godunov_include_diff_in_forcing)
            for (int lev = 0; lev <= finest_level; ++lev)
                MultiFab::Add(*vel_forces[lev], m_leveldata[lev]->divtau_o, 0, 0, AMREX_SPACEDIM, 0);

        if (nghost_force() > 0)
            fillpatch_force(m_cur_time, vel_forces, nghost_force());
    }

    // This will hold (1/rho) on faces
    Vector<MultiFab> inv_rho_x(finest_level+1);
    Vector<MultiFab> inv_rho_y(finest_level+1);
#if (AMREX_SPACEDIM == 3)
    Vector<MultiFab> inv_rho_z(finest_level+1);
#endif

    Vector<Array<MultiFab*,AMREX_SPACEDIM> > inv_rho(finest_level+1);
    for (int lev=0; lev <= finest_level; ++lev)
    {
        AMREX_D_TERM(inv_rho[lev][0] = &inv_rho_x[lev];,
                     inv_rho[lev][1] = &inv_rho_y[lev];,
                     inv_rho[lev][2] = &inv_rho_z[lev];);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        AMREX_D_TERM(
           inv_rho_x[lev].define(u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),Factory(lev));,
           inv_rho_y[lev].define(v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),Factory(lev));,
           inv_rho_z[lev].define(w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),Factory(lev)););

#ifdef AMREX_USE_EB
        EB_interp_CellCentroid_to_FaceCentroid (*density[lev], inv_rho[lev],
                                                0, 0, 1, geom[lev], get_density_bcrec());
#else
        amrex::average_cellcenter_to_face(inv_rho[lev], *density[lev], geom[lev]);
#endif

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            inv_rho[lev][idim]->invert(1.0, 0);
        }
    }

    compute_MAC_projected_velocities(vel, 
                                     AMREX_D_DECL(u_mac,v_mac,w_mac),
                                     AMREX_D_DECL(GetVecOfPtrs(inv_rho_x), GetVecOfPtrs(inv_rho_y),
                                                  GetVecOfPtrs(inv_rho_z)),
                                     vel_forces, time);

    // We now re-compute the velocity forcing terms including the pressure gradient,
    //    and compute the tracer forcing terms for the first time
    if (m_advection_type != "MOL")
    {
        compute_vel_forces(vel_forces, vel, density, tracer, tracer);

        if (m_godunov_include_diff_in_forcing)
            for (int lev = 0; lev <= finest_level; ++lev)
                MultiFab::Add(*vel_forces[lev], m_leveldata[lev]->divtau_o, 0, 0, AMREX_SPACEDIM, 0);

        if (nghost_force() > 0) 
            fillpatch_force(m_cur_time, vel_forces, nghost_force());

        // Note this is forcing for (rho s), not for s
        if (m_advect_tracer)
        {
            compute_tra_forces(tra_forces, get_density_old_const());
            if (m_godunov_include_diff_in_forcing)
                for (int lev = 0; lev <= finest_level; ++lev)
                    MultiFab::Add(*tra_forces[lev], m_leveldata[lev]->laps_o, 0, 0, m_ntrac, 0);
            if (nghost_force() > 0) 
                fillpatch_force(m_cur_time, tra_forces, nghost_force());
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (ngmac > 0) {
            AMREX_D_TERM(u_mac[lev]->FillBoundary(geom[lev].periodicity());,
                         v_mac[lev]->FillBoundary(geom[lev].periodicity());,
                         w_mac[lev]->FillBoundary(geom[lev].periodicity()););
        }

        MultiFab divu(vel[lev]->boxArray(),vel[lev]->DistributionMap(),1,2);
        divu.setVal(0.);
        Array<MultiFab const*, AMREX_SPACEDIM> u;
        u[0] = u_mac[lev];
        u[1] = v_mac[lev];
#if (AMREX_SPACEDIM == 3)
        u[2] = w_mac[lev];
#endif

#ifdef AMREX_USE_EB
        auto const& fact = EBFactory(lev);
        if (fact.isAllRegular())
            computeDivergence(divu,u,geom[lev]);
        else
            EB_computeDivergence(divu,u,geom[lev],true);
#else
        computeDivergence(divu,u,geom[lev]);
#endif
        divu.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        // Turn off tiling -- HACK HACK HACK
        for (MFIter mfi(*density[lev],false); mfi.isValid(); ++mfi)
        {

            Box const& bx = mfi.tilebox();
            compute_convective_term(bx, lev, mfi,
                                    conv_u[lev]->array(mfi),
                                    conv_r[lev]->array(mfi),
                                    (m_ntrac>0) ? conv_t[lev]->array(mfi) : Array4<Real>{},
                                    vel[lev]->const_array(mfi),
                                    density[lev]->array(mfi),
                                    (m_ntrac>0) ? tracer[lev]->const_array(mfi) : Array4<Real const>{},
                                    divu.const_array(mfi),
                                    AMREX_D_DECL(u_mac[lev]->const_array(mfi),
                                                 v_mac[lev]->const_array(mfi),
                                                 w_mac[lev]->const_array(mfi)),
                                    (!vel_forces.empty()) ? vel_forces[lev]->const_array(mfi)
                                                          : Array4<Real const>{},
                                    (!tra_forces.empty()) ? tra_forces[lev]->const_array(mfi)
                                                          : Array4<Real const>{});
        }
    }
}

void
incflo::compute_convective_term (Box const& bx, int lev, MFIter const& mfi,
                                 Array4<Real> const& dvdt, // velocity
                                 Array4<Real> const& drdt, // density
                                 Array4<Real> const& dtdt, // tracer
                                 Array4<Real const> const& vel,
                                 Array4<Real const> const& rho,
                                 Array4<Real const> const& tra,
                                 Array4<Real const> const& divu,
                                 AMREX_D_DECL(Array4<Real const> const& umac,
                                              Array4<Real const> const& vmac,
                                              Array4<Real const> const& wmac),
                                 Array4<Real const> const& fvel,
                                 Array4<Real const> const& ftra)
{
    Real l_dt = m_dt;

#ifdef AMREX_USE_EB
    auto const& fact = EBFactory(lev);
    EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
    Array4<EBCellFlag const> const& flag = flagfab.const_array();
    if (flagfab.getType(bx) == FabType::covered)
    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            AMREX_D_TERM(dvdt(i,j,k,0) = 0.0;,
                         dvdt(i,j,k,1) = 0.0;,
                         dvdt(i,j,k,2) = 0.0;);
            drdt(i,j,k) = 0.0;
        });
        if (m_advect_tracer) {
            amrex::ParallelFor(bx, m_ntrac, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dtdt(i,j,k,n) = 0.0;
            });
        }
        return;
    }

    bool regular = (flagfab.getType(amrex::grow(bx,2)) == FabType::regular);

    Array4<Real const> AMREX_D_DECL(fcx, fcy, fcz), ccc, vfrac, AMREX_D_DECL(apx, apy, apz);
    if (!regular) {
        AMREX_D_TERM(fcx = fact.getFaceCent()[0]->const_array(mfi);,
                     fcy = fact.getFaceCent()[1]->const_array(mfi);,
                     fcz = fact.getFaceCent()[2]->const_array(mfi););
        ccc   = fact.getCentroid().const_array(mfi);
        AMREX_D_TERM(apx = fact.getAreaFrac()[0]->const_array(mfi);,
                     apy = fact.getAreaFrac()[1]->const_array(mfi);,
                     apz = fact.getAreaFrac()[2]->const_array(mfi););
        vfrac = fact.getVolFrac().const_array(mfi);
    }
#endif

    Box rhotrac_box = amrex::grow(bx,2);
    if (m_advection_type != "MOL")  rhotrac_box.grow(1);
#ifdef AMREX_USE_EB
    if (!regular) rhotrac_box.grow(2);
#endif

    FArrayBox rhotracfab;
    Elixir eli_rt;
    Array4<Real> rhotrac;
    if (m_advect_tracer) {
        rhotracfab.resize(rhotrac_box, m_ntrac);
        eli_rt  = rhotracfab.elixir();
        rhotrac = rhotracfab.array();
        amrex::ParallelFor(rhotrac_box, m_ntrac,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            rhotrac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
        });
    }

    int nmaxcomp = AMREX_SPACEDIM;
    if (m_advect_tracer) nmaxcomp = std::max(nmaxcomp,m_ntrac);

    if (m_advection_type == "Godunov")
    {
        int n_tmp_fac;
        int n_tmp_grow;
#if (AMREX_SPACEDIM == 3)
        n_tmp_fac = 14;
#else
        n_tmp_fac = 10;
#endif

#ifdef AMREX_USE_EB
        n_tmp_grow = 4;
#else
        n_tmp_grow = 1;
#endif

        FArrayBox tmpfab(amrex::grow(bx,n_tmp_grow), nmaxcomp*n_tmp_fac+1);
        Elixir eli = tmpfab.elixir();

#ifdef AMREX_USE_EB
        Box gbx = bx;
        if (!regular)  
        {
            if (m_advection_type == "MOL") 
                gbx.grow(2);
            else if (m_advection_type == "Godunov") 
                gbx.grow(2);
            else 
                amrex::Abort("Dont know this advection_type");
        }
        // This one holds the convective term on a grown region so we can redistribute
        FArrayBox dUdt_tmpfab(gbx,nmaxcomp);
        Array4<Real> dUdt_tmp = dUdt_tmpfab.array();
        Elixir eli_du = dUdt_tmpfab.elixir();

        FArrayBox scratch_fab(gbx,3*nmaxcomp);
        Array4<Real> scratch = scratch_fab.array();
        Elixir eli_scratch = scratch_fab.elixir();

        if (!regular)
        {
            ebgodunov::compute_godunov_advection(gbx, AMREX_SPACEDIM,
                                                 dUdt_tmp, vel,
                                                 AMREX_D_DECL(umac, vmac, wmac), 
                                                 fvel, divu, l_dt, 
                                                 get_velocity_bcrec(),
                                                 get_velocity_bcrec_device_ptr(),
                                                 get_velocity_iconserv_device_ptr(),
                                                 tmpfab.dataPtr(), flag, 
                                                 AMREX_D_DECL(apx, apy, apz), vfrac,
                                                 AMREX_D_DECL(fcx, fcy, fcz), ccc, 
                                                 geom[lev], true); // is_velocity
            redistribution::redistribute_eb(bx, AMREX_SPACEDIM, dvdt, dUdt_tmp, vel, scratch,
                                            AMREX_D_DECL(umac, vmac, wmac), flag,
                                            AMREX_D_DECL(apx, apy, apz), vfrac,
                                            AMREX_D_DECL(fcx, fcy, fcz), ccc, geom[lev], l_dt, m_redistribution_type);

            if (!m_constant_density) {
                ebgodunov::compute_godunov_advection(gbx, 1,
                                                     dUdt_tmp, rho,
                                                     AMREX_D_DECL(umac, vmac, wmac), 
                                                     {}, divu, l_dt, 
                                                     get_density_bcrec(),
                                                     get_density_bcrec_device_ptr(),
                                                     get_density_iconserv_device_ptr(),
                                                     tmpfab.dataPtr(), flag,
                                                     AMREX_D_DECL(apx, apy, apz), vfrac,
                                                     AMREX_D_DECL(fcx, fcy, fcz), ccc, 
                                                     geom[lev]);
                redistribution::redistribute_eb(bx, 1, drdt, dUdt_tmp, rho, scratch,
                                                AMREX_D_DECL(umac, vmac, wmac), flag,
                                                AMREX_D_DECL(apx, apy, apz), vfrac,
                                                AMREX_D_DECL(fcx, fcy, fcz), ccc, geom[lev], l_dt, m_redistribution_type);
            }
            if (m_advect_tracer) {
                ebgodunov::compute_godunov_advection(gbx, m_ntrac,
                                                     dUdt_tmp, rhotrac,
                                                     AMREX_D_DECL(umac, vmac, wmac), 
                                                     ftra, divu, l_dt, 
                                                     get_tracer_bcrec(),
                                                     get_tracer_bcrec_device_ptr(),
                                                     get_tracer_iconserv_device_ptr(),
                                                     tmpfab.dataPtr(), flag,
                                                     AMREX_D_DECL(apx, apy, apz), vfrac,
                                                     AMREX_D_DECL(fcx, fcy, fcz), ccc, 
                                                     geom[lev]);
                redistribution::redistribute_eb(bx, m_ntrac, dtdt, dUdt_tmp, rhotrac, scratch,
                                                AMREX_D_DECL(umac, vmac, wmac), flag,
                                                AMREX_D_DECL(apx, apy, apz), vfrac,
                                                AMREX_D_DECL(fcx, fcy, fcz), ccc, geom[lev], l_dt, m_redistribution_type);
            }
            Gpu::streamSynchronize();
        }
        else
#endif
        {
            godunov::compute_godunov_advection(bx, AMREX_SPACEDIM,
                                               dvdt, vel,
                                               AMREX_D_DECL(umac, vmac, wmac), 
                                               fvel, divu, l_dt, 
                                               get_velocity_bcrec_device_ptr(),
                                               get_velocity_iconserv_device_ptr(),
                                               tmpfab.dataPtr(),m_godunov_ppm, 
                                               m_godunov_use_forces_in_trans, 
                                               geom[lev], true);
            if (!m_constant_density) {
                godunov::compute_godunov_advection(bx, 1,
                                                   drdt, rho,
                                                   AMREX_D_DECL(umac, vmac, wmac),
                                                   {}, divu, l_dt, 
                                                   get_density_bcrec_device_ptr(),
                                                   get_density_iconserv_device_ptr(),
                                                   tmpfab.dataPtr(),m_godunov_ppm,
                                                   m_godunov_use_forces_in_trans,
                                                   geom[lev]);
            }
            if (m_advect_tracer) {
                godunov::compute_godunov_advection(bx, m_ntrac,
                                                   dtdt, rhotrac,
                                                   AMREX_D_DECL(umac, vmac, wmac),
                                                   ftra, divu, l_dt, 
                                                   get_tracer_bcrec_device_ptr(),
                                                   get_tracer_iconserv_device_ptr(),
                                                   tmpfab.dataPtr(),m_godunov_ppm,
                                                   m_godunov_use_forces_in_trans,
                                                   geom[lev]);
            }
            Gpu::streamSynchronize();
        }
    }
    else if (m_advection_type == "MOL")
    {
        Box tmpbox = amrex::surroundingNodes(bx);
        int tmpcomp = nmaxcomp*3;
#ifdef AMREX_USE_EB
        Box gbx = bx;
        if (!regular) {
            gbx.grow(2);
            tmpbox.grow(3);
            tmpcomp += nmaxcomp;
        }
#endif
        FArrayBox tmpfab(tmpbox, tmpcomp);
        Elixir eli = tmpfab.elixir();

        AMREX_D_TERM(Array4<Real> fx = tmpfab.array(0);,
                     Array4<Real> fy = tmpfab.array(nmaxcomp);,
                     Array4<Real> fz = tmpfab.array(nmaxcomp*2););

#ifdef AMREX_USE_EB 
        if (!regular)
        {
            Array4<Real> scratch = tmpfab.array(0);
            Array4<Real> dUdt_tmp = tmpfab.array(nmaxcomp*3);

            // velocity
            mol::compute_convective_fluxes_eb(gbx, AMREX_SPACEDIM,
                                              AMREX_D_DECL(fx, fy, fz), vel, 
                                              AMREX_D_DECL(umac, vmac, wmac),
                                              get_velocity_bcrec().data(),
                                              get_velocity_bcrec_device_ptr(),
                                              flag, AMREX_D_DECL(fcx, fcy, fcz), ccc, geom[lev]);
            mol::compute_convective_rate_eb(gbx, AMREX_SPACEDIM, dUdt_tmp, AMREX_D_DECL(fx, fy, fz),
                                            flag, vfrac, AMREX_D_DECL(apx, apy, apz), geom[lev]);
            redistribution::redistribute_eb(bx, AMREX_SPACEDIM, dvdt, dUdt_tmp, vel, scratch,
                                            AMREX_D_DECL(umac, vmac, wmac), flag,
                                            AMREX_D_DECL(apx, apy, apz), vfrac,
                                            AMREX_D_DECL(fcx, fcy, fcz), ccc, geom[lev], l_dt, m_redistribution_type);

            // density
            if (!m_constant_density) {
                mol::compute_convective_fluxes_eb(gbx, 1,
                                                  AMREX_D_DECL(fx, fy, fz), rho, 
                                                  AMREX_D_DECL(umac, vmac, wmac),
                                                  get_density_bcrec().data(),
                                                  get_density_bcrec_device_ptr(),
                                                  flag, AMREX_D_DECL(fcx, fcy, fcz), ccc, geom[lev]);
                mol::compute_convective_rate_eb( gbx, 1, dUdt_tmp, AMREX_D_DECL(fx, fy, fz),
                                                flag, vfrac, AMREX_D_DECL(apx, apy, apz), geom[lev]);
                redistribution::redistribute_eb(bx, 1, drdt, dUdt_tmp, rho, scratch,
                                                AMREX_D_DECL(umac, vmac, wmac), flag,
                                                AMREX_D_DECL(apx, apy, apz), vfrac,
                                                AMREX_D_DECL(fcx, fcy, fcz), ccc, geom[lev], l_dt, m_redistribution_type);
            }

            if (m_advect_tracer) {
                mol::compute_convective_fluxes_eb(gbx, m_ntrac,
                                                  AMREX_D_DECL(fx, fy, fz), rhotrac, 
                                                  AMREX_D_DECL(umac, vmac, wmac),
                                                  get_tracer_bcrec().data(),
                                                  get_tracer_bcrec_device_ptr(),
                                                  flag, AMREX_D_DECL(fcx, fcy, fcz), ccc, geom[lev]);
                mol::compute_convective_rate_eb(gbx, m_ntrac, dUdt_tmp, AMREX_D_DECL(fx, fy, fz),
                                                flag, vfrac, AMREX_D_DECL(apx, apy, apz), geom[lev]);
                redistribution::redistribute_eb(bx, m_ntrac, dtdt, dUdt_tmp, rhotrac, scratch,
                                                AMREX_D_DECL(umac, vmac, wmac), flag,
                                                AMREX_D_DECL(apx, apy, apz), vfrac,
                                                AMREX_D_DECL(fcx, fcy, fcz), ccc, geom[lev], l_dt, m_redistribution_type);
            }
        }
        else
#endif
        {
            // const auto dxinv = geom[lev].InvCellSizeArray();

            // velocity
            mol::compute_convective_fluxes(bx, AMREX_SPACEDIM, AMREX_D_DECL(fx, fy, fz), vel,
                                           AMREX_D_DECL(umac, vmac, wmac),
                                           get_velocity_bcrec().data(),
                                           get_velocity_bcrec_device_ptr(),geom[lev]);

            // amrex_compute_divergence(bx,dvdt,AMREX_D_DECL(fx,fy,fz),dxinv);
            mol::compute_convective_rate(bx, AMREX_SPACEDIM, dvdt, AMREX_D_DECL(fx, fy, fz), geom[lev]);

            // density
            if (!m_constant_density) {
                mol::compute_convective_fluxes(bx, 1, AMREX_D_DECL(fx, fy, fz), rho,
                                               AMREX_D_DECL(umac, vmac, wmac),
                                               get_density_bcrec().data(),
                                               get_density_bcrec_device_ptr(),geom[lev]);
                mol::compute_convective_rate(bx, 1, drdt, AMREX_D_DECL(fx, fy, fz), geom[lev]);
            }

            // tracer
            if (m_advect_tracer) {
                mol::compute_convective_fluxes(bx, m_ntrac, AMREX_D_DECL(fx, fy, fz), rhotrac,
                                               AMREX_D_DECL(umac, vmac, wmac),
                                               get_tracer_bcrec().data(),
                                               get_tracer_bcrec_device_ptr(), geom[lev]);
                mol::compute_convective_rate(bx, m_ntrac, dtdt, AMREX_D_DECL(fx, fy, fz), geom[lev]);
            }
        }
    } else { 
        amrex::Abort("Dont know this advection type!");
    }
}

