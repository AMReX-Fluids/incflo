#include <incflo.H>
#include <prob_bc.H>
#include <hydro_godunov.H>
#include <hydro_mol.H>
#include <hydro_utils.H>
#include <AMReX_FillPatchUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#include <hydro_ebgodunov.H>
#include <hydro_ebmol.H>
#endif

using namespace amrex;

void incflo::init_advection ()
{
    // Use convective differencing for velocity
    if (m_advect_momentum) {
        m_iconserv_velocity.resize(  AMREX_SPACEDIM, 1);
        m_iconserv_velocity_d.resize(AMREX_SPACEDIM, 1);
    } else {
        m_iconserv_velocity.resize(  AMREX_SPACEDIM, 0);
        m_iconserv_velocity_d.resize(AMREX_SPACEDIM, 0);
    }

    // Density is always updated conservatively
    m_iconserv_density.resize(1, 1);
    m_iconserv_density_d.resize(1, 1);

    // We update (rho * tracer), not tracer itself, hence we update conservatively
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
    bool fluxes_are_area_weighted = false;
    bool knownFaceStates          = false; // HydroUtils always recompute face states

#ifdef AMREX_USE_EB
    amrex::Print() << "REDISTRIBUTION TYPE " << m_redistribution_type << std::endl;
#endif

    // Make one flux MF at each level to hold all the fluxes (velocity, density, tracers)
    int n_flux_comp = AMREX_SPACEDIM;
    if (!m_constant_density) n_flux_comp += 1;
    if ( m_advect_tracer)    n_flux_comp += m_ntrac;

    // This will hold state on faces
    Vector<MultiFab> face_x(finest_level+1);
    Vector<MultiFab> face_y(finest_level+1);
#if (AMREX_SPACEDIM == 3)
    Vector<MultiFab> face_z(finest_level+1);
#endif

    // This will hold fluxes on faces
    Vector<MultiFab> flux_x(finest_level+1);
    Vector<MultiFab> flux_y(finest_level+1);
#if (AMREX_SPACEDIM == 3)
    Vector<MultiFab> flux_z(finest_level+1);
#endif

    Vector<MultiFab> divu(finest_level+1);
    Vector<MultiFab> rhovel(finest_level+1);
    Vector<MultiFab> rhotrac(finest_level+1);

    Vector<Array<MultiFab*,AMREX_SPACEDIM> > fluxes(finest_level+1);
    Vector<Array<MultiFab*,AMREX_SPACEDIM> >  faces(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(
           face_x[lev].define(u_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),u_mac[lev]->Factory());,
           face_y[lev].define(v_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),v_mac[lev]->Factory());,
           face_z[lev].define(w_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),w_mac[lev]->Factory()););
        AMREX_D_TERM(
           flux_x[lev].define(u_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),u_mac[lev]->Factory());,
           flux_y[lev].define(v_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),v_mac[lev]->Factory());,
           flux_z[lev].define(w_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),w_mac[lev]->Factory()););

        divu[lev].define(vel[lev]->boxArray(),dmap[lev],1,4,MFInfo(),vel[lev]->Factory());
        if (m_advect_momentum)
            rhovel[lev].define(vel[lev]->boxArray(),dmap[lev],AMREX_SPACEDIM,
                               vel[lev]->nGrow(),MFInfo(),vel[lev]->Factory());
        if (m_advect_tracer && m_ntrac > 0)
            rhotrac[lev].define(tracer[lev]->boxArray(),dmap[lev],tracer[lev]->nComp(),
                                tracer[lev]->nGrow(),MFInfo(),tracer[lev]->Factory());

        AMREX_D_TERM(faces[lev][0] = &face_x[lev];,
                     faces[lev][1] = &face_y[lev];,
                     faces[lev][2] = &face_z[lev];);

        AMREX_D_TERM(fluxes[lev][0] = &flux_x[lev];,
                     fluxes[lev][1] = &flux_y[lev];,
                     fluxes[lev][2] = &flux_z[lev];);
    }

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
        if (nghost_mac() > 0)
        {
            // FillPatch umac.
            Array<MultiFab*, AMREX_SPACEDIM> u_fine;
            AMREX_D_TERM(u_fine[0] = u_mac[lev  ];,
                         u_fine[1] = v_mac[lev  ];,
                         u_fine[2] = w_mac[lev  ];);

            PhysBCFunct<GpuBndryFuncFab<IncfloVelFill>>
                fine_bndry_func_x(geom[lev],
                                  m_bcrec_velocity,
                                  IncfloVelFill{m_probtype, m_bc_velocity});
            PhysBCFunct<GpuBndryFuncFab<IncfloVelFill>>
                fine_bndry_func_y(geom[lev],
                                  m_bcrec_velocity,
                                  IncfloVelFill{m_probtype, m_bc_velocity});
#if (AMREX_SPACEDIM == 3)
            PhysBCFunct<GpuBndryFuncFab<IncfloVelFill>>
                fine_bndry_func_z(geom[lev],
                                  m_bcrec_velocity,
                                  IncfloVelFill{m_probtype, m_bc_velocity});
#endif
            Array<PhysBCFunct<GpuBndryFuncFab<IncfloVelFill>>,AMREX_SPACEDIM>
                fbndyFuncArr = {AMREX_D_DECL(fine_bndry_func_x,fine_bndry_func_y,fine_bndry_func_z)};

            if (lev == 0)
            {
                //
                // BDS needs umac on physical boundaries.
                // Godunov handles physical boundaries interally, but needs periodic ghosts filled.
                // MOL doesn't need any umac ghost cells, so it doesn't get here.
                //
                Real fake_time = 0.;

                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    amrex::FillPatchSingleLevel(*u_fine[idim], IntVect(nghost_mac()), fake_time,
                                                {u_fine[idim]}, {fake_time},
                                                0, 0, 1, geom[lev],
                                                fbndyFuncArr[idim], idim);
                }
            }
            else // lev > 0
            {
                IntVect rr  = geom[lev].Domain().size() / geom[lev-1].Domain().size();
                Array<MultiFab*, AMREX_SPACEDIM> u_crse;
                AMREX_D_TERM(u_crse[0] = u_mac[lev-1];,
                             u_crse[1] = v_mac[lev-1];,
                             u_crse[2] = w_mac[lev-1];);


                // Divergence preserving interp
                Interpolater* mapper = &face_divfree_interp;

                const Array<Vector<BCRec>,AMREX_SPACEDIM> bcrecArr = {AMREX_D_DECL(m_bcrec_velocity,
                                                                                   m_bcrec_velocity,
                                                                                   m_bcrec_velocity)};

                PhysBCFunct<GpuBndryFuncFab<IncfloVelFill>>
                    crse_bndry_func(geom[lev-1],
                                    m_bcrec_velocity,
                                    IncfloVelFill{m_probtype, m_bc_velocity});
                Array<PhysBCFunct<GpuBndryFuncFab<IncfloVelFill>>,AMREX_SPACEDIM>
                    cbndyFuncArr = {AMREX_D_DECL(crse_bndry_func,crse_bndry_func,crse_bndry_func)};

                // Use piecewise constant interpolation in time, so create dummy variable for time
                Real fake_time = 0.;
                Array<int, AMREX_SPACEDIM> idx = {AMREX_D_DECL(0,1,2)};
                FillPatchTwoLevels(u_fine, IntVect(nghost_mac()), fake_time,
                                   {u_crse}, {fake_time},
                                   {u_fine}, {fake_time},
                                   0, 0, 1,
                                   geom[lev-1], geom[lev],
                                   cbndyFuncArr, idx, fbndyFuncArr, idx,
                                   rr, mapper, bcrecArr, idx);
            }
        } // end umac fill

        divu[lev].setVal(0.);
        Array<MultiFab const*, AMREX_SPACEDIM> u;
        AMREX_D_TERM(u[0] = u_mac[lev];,
                     u[1] = v_mac[lev];,
                     u[2] = w_mac[lev];);

        // // CEG fixme
        // VisMF::Write(*u_mac[0],"umac");
        // VisMF::Write(*v_mac[0],"vmac");
        // VisMF::Write(*vel[0],"vpred");
        // VisMF::Write(*density[0],"rpred2");

#ifdef AMREX_USE_EB
        const auto& ebfact_old = OldEBFactory(lev);
        const auto& ebfact_new =    EBFactory(lev);
        const auto& ebfact =    EBFactory(lev, time);

        if (!ebfact_old.isAllRegular()) {
            // NOTE this divu is not relevant for MSRD
            // It's used for non-conservative adjustment (but MSRD is always conservative)
            // and in EBGod in 3D for corner-coupling (but flow through EB doesn't use any
            // transverse or CC (aka d/dt) terms)
            if (m_eb_flow.enabled) {
                amrex::EB_computeDivergence(divu[lev],u,geom[lev],true,*get_velocity_eb()[lev]);
            } else {
                amrex::EB_computeDivergence(divu[lev],u,geom[lev],true);
            }
        }
        else
#endif
        {
            amrex::computeDivergence(divu[lev],u,geom[lev]);
        }

        divu[lev].FillBoundary(geom[lev].periodicity());

        // ************************************************************************
        // Compute advective fluxes
        // ************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*density[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();

            Array4<Real const> const& divu_arr = divu[lev].const_array(mfi);

            // ************************************************************************
            // Velocity
            // ************************************************************************
            // Not sure if this temporary for rho*forces is the best option...
            FArrayBox rhovel_f;
            if ( m_advect_momentum )
            {
                // create rho*U

                // Note we must actually grow the tilebox, not use growntilebox, because
                // we want to use this immediately below and we need all the "ghost cells" of
                // the tiled region
                Box const& bxg = amrex::grow(bx,vel[lev]->nGrow());

                Array4<Real const> U       =     vel[lev]->const_array(mfi);
                Array4<Real const> rho     = density[lev]->const_array(mfi);
                Array4<Real      > rho_vel =  rhovel[lev].array(mfi);

                amrex::ParallelFor(bxg, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rho_vel(i,j,k,n) = rho(i,j,k) * U(i,j,k,n);
                });

                if (!vel_forces.empty())
                {
                    Box const& fbx = amrex::grow(bx,nghost_force());
                    rhovel_f.resize(fbx, AMREX_SPACEDIM, The_Async_Arena());
                    Array4<Real const> vf        = vel_forces[lev]->const_array(mfi);
                    Array4<Real      > rho_vel_f =  rhovel_f.array();

                    amrex::ParallelFor(fbx, AMREX_SPACEDIM,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        rho_vel_f(i,j,k,n) = rho(i,j,k) * vf(i,j,k,n);
                    });
                }
            }

            int face_comp = 0;
            int ncomp = AMREX_SPACEDIM;
            bool is_velocity = true;
            HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi,
                                                     (m_advect_momentum) ? rhovel[lev].array(mfi) : vel[lev]->const_array(mfi),
                                                     AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                                  flux_y[lev].array(mfi,face_comp),
                                                                  flux_z[lev].array(mfi,face_comp)),
                                                     AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                                  face_y[lev].array(mfi,face_comp),
                                                                  face_z[lev].array(mfi,face_comp)),
                                                     knownFaceStates,
                                                     AMREX_D_DECL(u_mac[lev]->const_array(mfi),
                                                                  v_mac[lev]->const_array(mfi),
                                                                  w_mac[lev]->const_array(mfi)),
                                                     divu_arr,
                                                     (vel_forces.empty())
                                                         ? Array4<Real const>{}
                                                         : m_advect_momentum
                                                             ? rhovel_f.const_array()
                                                             : vel_forces[lev]->const_array(mfi),
                                                     geom[lev], m_dt,
                                                     get_velocity_bcrec(),
                                                     get_velocity_bcrec_device_ptr(),
                                                     get_velocity_iconserv_device_ptr(),
#ifdef AMREX_USE_EB
                                                     ebfact,
                                                     m_eb_flow.enabled ? get_velocity_eb()[lev]->const_array(mfi) : Array4<Real const>{},
#endif
                                                     m_godunov_ppm, m_godunov_use_forces_in_trans,
                                                     is_velocity, fluxes_are_area_weighted,
                                                     m_advection_type);


            // ************************************************************************
            // Density
            // ************************************************************************
            if (!m_constant_density)
            {
                face_comp = AMREX_SPACEDIM;
                ncomp = 1;
                is_velocity = false;
                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi,
                                                         density[lev]->const_array(mfi),
                                                         AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                                      flux_y[lev].array(mfi,face_comp),
                                                                      flux_z[lev].array(mfi,face_comp)),
                                                         AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                                      face_y[lev].array(mfi,face_comp),
                                                                      face_z[lev].array(mfi,face_comp)),
                                                         knownFaceStates,
                                                         AMREX_D_DECL(u_mac[lev]->const_array(mfi),
                                                                      v_mac[lev]->const_array(mfi),
                                                                      w_mac[lev]->const_array(mfi)),
                                                         divu_arr, Array4<Real const>{},
                                                         geom[lev], m_dt,
                                                         get_density_bcrec(),
                                                         get_density_bcrec_device_ptr(),
                                                         get_density_iconserv_device_ptr(),
#ifdef AMREX_USE_EB
                                                         ebfact,
                                                         m_eb_flow.enabled ? get_density_eb()[lev]->const_array(mfi) : Array4<Real const>{},
#endif
                                                         m_godunov_ppm, m_godunov_use_forces_in_trans,
                                                         is_velocity, fluxes_are_area_weighted,
                                                         m_advection_type);
            }

            // ************************************************************************
            // (Rho*Tracer)
            // ************************************************************************
            // Make a FAB holding (rho * tracer) that is the same size as the original tracer FAB
            if (m_advect_tracer && (m_ntrac>0)) {

                // Note we must actually grow the tilebox, not use growntilebox, because
                // we want to use this immediately below and we need all the "ghost cells" of
                // the tiled region
                Box const& bxg = amrex::grow(bx,tracer[lev]->nGrow());

                Array4<Real const> tra     =  tracer[lev]->const_array(mfi);
                Array4<Real const> rho     = density[lev]->const_array(mfi);
                Array4<Real      > ro_trac = rhotrac[lev].array(mfi);

                amrex::ParallelFor(bxg, m_ntrac,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    ro_trac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
                });

                if (m_constant_density)
                   face_comp = AMREX_SPACEDIM;
                else
                   face_comp = AMREX_SPACEDIM+1;
                ncomp = m_ntrac;
                is_velocity = false;

                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi, ro_trac,
                                          AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                       flux_y[lev].array(mfi,face_comp),
                                                       flux_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                       face_y[lev].array(mfi,face_comp),
                                                       face_z[lev].array(mfi,face_comp)),
                                          knownFaceStates,
                                          AMREX_D_DECL(u_mac[lev]->const_array(mfi),
                                                       v_mac[lev]->const_array(mfi),
                                                       w_mac[lev]->const_array(mfi)),
                                          divu_arr,
                                          (!tra_forces.empty()) ? tra_forces[lev]->const_array(mfi) : Array4<Real const>{},
                                          geom[lev], m_dt,
                                          get_tracer_bcrec(),
                                          get_tracer_bcrec_device_ptr(),
                                          get_tracer_iconserv_device_ptr(),
#ifdef AMREX_USE_EB
                                          ebfact,
                                          m_eb_flow.enabled ? get_tracer_eb()[lev]->const_array(mfi) : Array4<Real const>{},
#endif
                                          m_godunov_ppm, m_godunov_use_forces_in_trans,
                                          is_velocity, fluxes_are_area_weighted,
                                          m_advection_type);
            }
        } // mfi
    } // lev

    // In order to enforce conservation across coarse-fine boundaries we must be sure to average down the fluxes
    //    before we use them.  Note we also need to average down the face states if we are going to do
    //    convective differencing
    for (int lev = finest_level; lev > 0; --lev)
    {
        IntVect rr  = geom[lev].Domain().size() / geom[lev-1].Domain().size();
#ifdef AMREX_USE_EB
        EB_average_down_faces(GetArrOfConstPtrs( faces[lev]),  faces[lev-1], rr, geom[lev-1]);
        EB_average_down_faces(GetArrOfConstPtrs(fluxes[lev]), fluxes[lev-1], rr, geom[lev-1]);
#else
        average_down_faces(GetArrOfConstPtrs( faces[lev]),  faces[lev-1], rr, geom[lev-1]);
        average_down_faces(GetArrOfConstPtrs(fluxes[lev]), fluxes[lev-1], rr, geom[lev-1]);
#endif
    }

    for (int lev = 0; lev <= finest_level; ++lev)
    {
#ifdef AMREX_USE_EB
        // FIXME - need to think about how to make this work for not Moving EB
        const EBFArrayBoxFactory* ebfact_new = &EBFactory(lev, m_cur_time+m_dt);
        const EBFArrayBoxFactory* ebfact     = &EBFactory(lev, time);
        auto const& vfrac = ebfact->getVolFrac();

        // //fixme
        // VisMF::Write(vfrac, "vfm");
        // VisMF::Write(ebfact_new->getVolFrac(), "vfn");
        // VisMF::Write(*get_velocity_eb()[0], "veb");

        // MSRD updates are really associated to the EB at both times...
        // Don't think it actaully matters what time this factory is at though.
        int tmp_ng = 3;
        MultiFab dvdt_tmp(vel[lev]->boxArray(),dmap[lev],AMREX_SPACEDIM,tmp_ng,MFInfo(),Factory(lev));
        MultiFab drdt_tmp(vel[lev]->boxArray(),dmap[lev],1             ,tmp_ng,MFInfo(),Factory(lev));
        MultiFab dtdt_tmp(vel[lev]->boxArray(),dmap[lev],m_ntrac       ,tmp_ng,MFInfo(),Factory(lev));

        // Must initialize to zero because not all values may be set, e.g. outside the domain.
        dvdt_tmp.setVal(0.);
        drdt_tmp.setVal(0.);
        dtdt_tmp.setVal(0.);
#endif

        Real mult = -1.0;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*conv_u[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();

            if (m_verbose > 1) { Print()<<"Vel advection term..."<<std::endl; }

            int flux_comp = 0;
            int  num_comp = AMREX_SPACEDIM;
#ifdef AMREX_USE_EB
            FArrayBox rhovel_eb;
            if (m_advect_momentum && m_eb_flow.enabled)
            {
                // FIXME - not sure if rhovel_eb needs same num grow cells as vel or if
                // could make do with tmp_ng
                Box const& bxg = amrex::grow(bx,vel[lev]->nGrow());
                rhovel_eb.resize(bxg, AMREX_SPACEDIM, The_Async_Arena());

                Array4<Real const> U       = get_velocity_eb()[lev]->const_array(mfi);
                Array4<Real const> rho     =  get_density_eb()[lev]->const_array(mfi);
                Array4<Real      > rho_vel =  rhovel_eb.array();

                amrex::ParallelFor(bxg, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rho_vel(i,j,k,n) = rho(i,j,k) * U(i,j,k,n);
                });
            }
            EBCellFlagFab const& flagfab = ebfact->getMultiEBCellFlagFab()[mfi];
            auto const& update_arr  = dvdt_tmp.array(mfi);
            if (flagfab.getType(bx) != FabType::covered)
                HydroUtils::EB_ComputeDivergence(bx, update_arr,
                                                 AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                              flux_y[lev].const_array(mfi,flux_comp),
                                                              flux_z[lev].const_array(mfi,flux_comp)),
                                                 vfrac.const_array(mfi), num_comp, geom[lev],
                                                 mult, fluxes_are_area_weighted,
#ifdef AMREX_USE_MOVING_EB
// Don't flux through moving EB. EB velocity in computing edgestates/fluxes only serves to prohibit using
// any d/dt terms (e.g transverse)
                                                 // For corrector, we include the contribution from EB
                                                 (time==m_cur_time) ? Array4<Real const>{} : get_velocity_eb()[lev]->const_array(mfi),
                                                 (time==m_cur_time) ? Array4<Real const>{}
                                                 : ( m_advect_momentum  ? rhovel_eb.const_array() : get_velocity_eb()[lev]->const_array(mfi)),
#else
                                                 m_eb_flow.enabled ?
                                                    get_velocity_eb()[lev]->const_array(mfi) : Array4<Real const>{},
                                                 m_eb_flow.enabled ?
                                                 ( m_advect_momentum  ? rhovel_eb.const_array() : get_velocity_eb()[lev]->const_array(mfi))
                                                 : Array4<Real const>{},
#endif
                                                 flagfab.const_array(),
                                                 (flagfab.getType(bx) != FabType::regular) ?
                                                    ebfact->getBndryArea().const_array(mfi) : Array4<Real const>{},
                                                 (flagfab.getType(bx) != FabType::regular) ?
                                                    ebfact->getBndryNormal().const_array(mfi) : Array4<Real const>{});
#else
            auto const& update_arr  = conv_u[lev]->array(mfi);
            HydroUtils::ComputeDivergence(bx, update_arr,
                                          AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                       flux_y[lev].const_array(mfi,flux_comp),
                                                       flux_z[lev].const_array(mfi,flux_comp)),
                                          num_comp, geom[lev],
                                          mult, fluxes_are_area_weighted);
#endif

            if (!m_advect_momentum)
            {
                // For convective, we define u dot grad u = div (u u) - u div(u)
                HydroUtils::ComputeConvectiveTerm(bx, num_comp, mfi,
                                                  vel[lev]->array(mfi,0),
                                                  AMREX_D_DECL(face_x[lev].array(mfi),
                                                               face_y[lev].array(mfi),
                                                               face_z[lev].array(mfi)),
                                                  divu[lev].array(mfi),
                                                  update_arr,
                                                  get_velocity_iconserv_device_ptr(),
#ifdef AMREX_USE_EB
                                                  *ebfact,
#endif
                                                  m_advection_type);
            }
        } // end mfi

        // Note: density is always updated conservatively -- we do not provide an option for
        //       updating density convectively
        if (!m_constant_density)
        {
          int flux_comp = AMREX_SPACEDIM;

          for (MFIter mfi(*conv_r[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
            if (m_verbose > 1) { Print()<<"Density advection term..."<<std::endl; }

            EBCellFlagFab const& flagfab = ebfact->getMultiEBCellFlagFab()[mfi];
            if (flagfab.getType(bx) != FabType::covered)
                HydroUtils::EB_ComputeDivergence(bx, drdt_tmp.array(mfi),
                                                 AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                              flux_y[lev].const_array(mfi,flux_comp),
                                                              flux_z[lev].const_array(mfi,flux_comp)),
                                                 vfrac.const_array(mfi), 1, geom[lev], mult,
                                                 fluxes_are_area_weighted,
#ifdef AMREX_USE_MOVING_EB
// Don't flux through moving EB. EB velocity in computing edgestates/fluxes only serves to prohibit using
// any d/dt terms (e.g transverse)
                                                 (time==m_cur_time) ? Array4<Real const>{} : get_velocity_eb()[lev]->const_array(mfi),
                                                 (time==m_cur_time) ? Array4<Real const>{} : get_density_eb()[lev]->const_array(mfi),
#else
                                                 m_eb_flow.enabled ?
                                                    get_velocity_eb()[lev]->const_array(mfi) : Array4<Real const>{},
                                                 m_eb_flow.enabled ?
                                                    get_density_eb()[lev]->const_array(mfi) : Array4<Real const>{},
#endif
                                                 flagfab.const_array(),
                                                 (flagfab.getType(bx) != FabType::regular) ?
                                                    ebfact->getBndryArea().const_array(mfi) : Array4<Real const>{},
                                                 (flagfab.getType(bx) != FabType::regular) ?
                                                    ebfact->getBndryNormal().const_array(mfi) : Array4<Real const>{});
#else
            HydroUtils::ComputeDivergence(bx, conv_r[lev]->array(mfi),
                                          AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                       flux_y[lev].const_array(mfi,flux_comp),
                                                       flux_z[lev].const_array(mfi,flux_comp)),
                                          1, geom[lev], mult,
                                          fluxes_are_area_weighted);
#endif
          } // mfi
        } // not constant density

        // Note: (rho*trac) is always updated conservatively -- we do not provide an option for
        //       updating (rho*trac) convectively
        if (m_advect_tracer && m_ntrac > 0)
        {
          int flux_comp = (m_constant_density) ? AMREX_SPACEDIM : AMREX_SPACEDIM+1;

          for (MFIter mfi(*conv_t[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
            EBCellFlagFab const& flagfab = ebfact->getMultiEBCellFlagFab()[mfi];
            auto const& update_arr  = dtdt_tmp.array(mfi);
            if (flagfab.getType(bx) != FabType::covered)
                HydroUtils::EB_ComputeDivergence(bx, update_arr,
                                                 AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                              flux_y[lev].const_array(mfi,flux_comp),
                                                              flux_z[lev].const_array(mfi,flux_comp)),
                                                 vfrac.const_array(mfi), m_ntrac, geom[lev], mult,
                                                 fluxes_are_area_weighted,

#ifdef AMREX_USE_MOVING_EB
// Don't flux through moving EB. EB velocity in computing edgestates/fluxes only serves to prohibit using
// any d/dt terms (e.g transverse)
                                                 // this didn't help with corrector
                                                 (time==m_cur_time) ? Array4<Real const>{} : get_velocity_eb()[lev]->const_array(mfi),
                                                 (time==m_cur_time) ? Array4<Real const>{} : get_tracer_eb()[lev]->const_array(mfi),
#else
                                                 m_eb_flow.enabled ?
                                                    get_velocity_eb()[lev]->const_array(mfi) : Array4<Real const>{},
                                                 m_eb_flow.enabled ?
                                                    get_tracer_eb()[lev]->const_array(mfi) : Array4<Real const>{},
#endif
                                                 flagfab.const_array(),
                                                 (flagfab.getType(bx) != FabType::regular) ?
                                                    ebfact->getBndryArea().const_array(mfi) : Array4<Real const>{},
                                                 (flagfab.getType(bx) != FabType::regular) ?
                                                    ebfact->getBndryNormal().const_array(mfi) : Array4<Real const>{});
#else
                auto const& update_arr  = conv_t[lev]->array(mfi);
                HydroUtils::ComputeDivergence(bx, update_arr,
                                              AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                           flux_y[lev].const_array(mfi,flux_comp),
                                                           flux_z[lev].const_array(mfi,flux_comp)),
                                              m_ntrac, geom[lev], mult,
                                              fluxes_are_area_weighted);
#endif

          } // mfi
        } // advect tracer

        // //fixme
        // VisMF::Write(dvdt_tmp,"vtmp");
        // VisMF::Write(drdt_tmp,"rtmp");

#ifdef AMREX_USE_EB
        // We only filled these on the valid cells so we fill same-level interior ghost cells here.
        // (We don't need values outside the domain or at a coarser level so we can call just FillBoundary)
        dvdt_tmp.FillBoundary(geom[lev].periodicity());
        drdt_tmp.FillBoundary(geom[lev].periodicity());
        dtdt_tmp.FillBoundary(geom[lev].periodicity());

//fixme
        // VisMF::Write(dvdt_tmp,"vtmp");
        // VisMF::Write(drdt_tmp,"rtmp");

        for (MFIter mfi(*density[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            if ( time == m_cur_time )
            {
            redistribute_convective_term (bx, mfi,
                                          (m_advect_momentum) ? rhovel[lev].const_array(mfi) : vel[lev]->const_array(mfi),
                                          density[lev]->const_array(mfi),
                                          (m_advect_tracer && (m_ntrac>0)) ? rhotrac[lev].const_array(mfi) : Array4<Real const>{},
                                          dvdt_tmp.array(mfi),
                                          drdt_tmp.array(mfi),
                                          (m_advect_tracer && (m_ntrac>0)) ? dtdt_tmp.array(mfi) : Array4<Real>{},
                                          conv_u[lev]->array(mfi),
                                          conv_r[lev]->array(mfi),
                                          (m_advect_tracer && (m_ntrac>0)) ? conv_t[lev]->array(mfi) : Array4<Real>{},
                                          m_redistribution_type, m_constant_density, m_advect_tracer, m_ntrac,
                                              ebfact, ebfact_new,
                                              m_eb_flow.enabled ? get_velocity_eb()[lev]->const_array(mfi) : Array4<Real const>{},
                                              geom[lev], m_dt);
            }
            else
            {
                // Regular SRD
                // When we turn SRD slopes on, we need to think more about what we really want to be
                // redistributing here
                redistribute_convective_term (bx, mfi,
                                              (m_advect_momentum) ? rhovel[lev].const_array(mfi) : vel[lev]->const_array(mfi),
                                              density[lev]->const_array(mfi),
                                              (m_advect_tracer && (m_ntrac>0)) ? rhotrac[lev].const_array(mfi) : Array4<Real const>{},
                                              dvdt_tmp.array(mfi),
                                              drdt_tmp.array(mfi),
                                              (m_advect_tracer && (m_ntrac>0)) ? dtdt_tmp.array(mfi) : Array4<Real>{},
                                              conv_u[lev]->array(mfi),
                                              conv_r[lev]->array(mfi),
                                              (m_advect_tracer && (m_ntrac>0)) ? conv_t[lev]->array(mfi) : Array4<Real>{},
                                              m_redistribution_type, m_constant_density, m_advect_tracer, m_ntrac,
                                              ebfact, ebfact, // FIXME, for now need these the same
                                              // m_eb_flow.enabled ? get_velocity_eb()[lev]->const_array(mfi) : Array4<Real const>{},
                                              Array4<Real const>{}, // No moving EB here, need a single switch...
                                              geom[lev], m_dt);
            }
       } // mfi
#endif
    } // lev
}
