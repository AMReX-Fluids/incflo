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

    // Advect scalars conservatively?
    m_iconserv_tracer.resize(m_ntrac, 1);
    ParmParse pp("incflo");
    pp.queryarr("trac_is_conservative", m_iconserv_tracer, 0, m_ntrac );
    m_iconserv_tracer_d.resize(m_ntrac);
    // copy
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy
#else
    std::memcpy
#endif
        (m_iconserv_tracer_d.data(), m_iconserv_tracer.data(), sizeof(int)*m_ntrac);

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
                                 Real /*time*/)
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

    bool any_conserv_trac = false;
    for (auto& i : m_iconserv_tracer){
        if ( i == 1 ) {
            any_conserv_trac = true;
            break;
        }
    }
    bool any_convective_trac = false;
    for (auto& i : m_iconserv_tracer){
        if ( i == 0 ) {
            any_convective_trac = true;
            break;
        }
    }

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

        if (m_advect_momentum) {
            rhovel[lev].define(vel[lev]->boxArray(),dmap[lev],AMREX_SPACEDIM,
                               vel[lev]->nGrow(),MFInfo(),vel[lev]->Factory());
        }
        if (m_advect_tracer && m_ntrac > 0 && any_conserv_trac) {
            rhotrac[lev].define(tracer[lev]->boxArray(),dmap[lev],tracer[lev]->nComp(),
                                tracer[lev]->nGrow(),MFInfo(),tracer[lev]->Factory());
        }

        AMREX_D_TERM(faces[lev][0] = &face_x[lev];,
                     faces[lev][1] = &face_y[lev];,
                     faces[lev][2] = &face_z[lev];);

        AMREX_D_TERM(fluxes[lev][0] = &flux_x[lev];,
                     fluxes[lev][1] = &flux_y[lev];,
                     fluxes[lev][2] = &flux_z[lev];);
    }

    // We now re-compute the velocity forcing terms including the pressure gradient,
    //    and compute the tracer forcing terms for the first time
    if (m_advection_type != "MOL") {

        compute_vel_forces(vel_forces, vel, density, tracer, tracer);

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

        if (nghost_force() > 0)
            fillpatch_force(m_cur_time, vel_forces, nghost_force());

        // Note that for conservative tracers, this is forcing for (rho s)
        // and for non-conservative, this is forcing for s
        if (m_advect_tracer)
        {
            compute_tra_forces(tra_forces, get_density_old_const());
            if (m_godunov_include_diff_in_forcing)
                for (int lev = 0; lev <= finest_level; ++lev)
                    MultiFab::Add(*tra_forces[lev], m_leveldata[lev]->laps_o, 0, 0, m_ntrac, 0);
            if (nghost_force() > 0)
                fillpatch_force(m_cur_time, tra_forces, nghost_force());
        }

    } // end m_advection_type

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
                // Godunov handles physical boundaries internally, but needs periodic ghosts filled.
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

#ifdef AMREX_USE_EB
        const auto& ebfact = EBFactory(lev);

        if (!ebfact.isAllRegular()) {
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

        // *************************************************************************************
        // Define domain boundary conditions at half-time to be used for fluxes if using Godunov
        // *************************************************************************************
        //
        MultiFab vel_nph, rho_nph, trac_nph;
        if (m_advection_type != "MOL") {
            vel_nph.define(vel[lev]->boxArray(),vel[lev]->DistributionMap(),AMREX_SPACEDIM,1);
            vel_nph.setVal(0.);
            fillphysbc_velocity(lev, m_cur_time+0.5*m_dt, vel_nph, 1);

            if ( !m_constant_density || m_advect_momentum ||
                (m_advect_tracer && m_ntrac > 0) )
            {
                rho_nph.define(density[lev]->boxArray(),density[lev]->DistributionMap(),1,1);
                rho_nph.setVal(0.);
                fillphysbc_density(lev, m_cur_time+0.5*m_dt, rho_nph, 1);
            }

            if ( m_advect_momentum ) {
                for (int n = 0; n < AMREX_SPACEDIM; n++) {
                    Multiply(vel_nph, rho_nph, 0, n, 1, 1);
                }
            }

            if (m_advect_tracer && (m_ntrac>0)) {
                trac_nph.define(tracer[lev]->boxArray(),tracer[lev]->DistributionMap(),m_ntrac,1);
                trac_nph.setVal(0.);
                fillphysbc_tracer(lev, m_cur_time+0.5*m_dt, trac_nph, 1);
                auto const* iconserv = get_tracer_iconserv_device_ptr();
                for (int n = 0; n < m_ntrac; n++) {
                    if ( iconserv[n] ){
                        Multiply(trac_nph, rho_nph, 0, n, 1, 1);
                    }
                }
            }
        }

        // ************************************************************************
        // Define mixed boundary conditions if relevant
        // ************************************************************************
        //
        // Create BC MF first (to hold bc's that vary in space along the face)
        //
        std::unique_ptr<iMultiFab> velBC_MF;
        if (m_has_mixedBC) {
            velBC_MF = make_BC_MF(lev, m_bcrec_velocity_d, "velocity");
        }
        std::unique_ptr<iMultiFab> densBC_MF;
        if (m_has_mixedBC) {
            densBC_MF = make_BC_MF(lev, m_bcrec_density_d, "density");
        }
        std::unique_ptr<iMultiFab> tracBC_MF;
        if (m_advect_tracer  && (m_ntrac>0)) {
            if (m_has_mixedBC) {
                tracBC_MF = make_BC_MF(lev, m_bcrec_tracer_d, "tracer");
            }
        }

        // ************************************************************************
        // Compute advective fluxes
        // ************************************************************************
        //
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

                ParallelFor(bxg, AMREX_SPACEDIM,
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

                    ParallelFor(fbx, AMREX_SPACEDIM,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        rho_vel_f(i,j,k,n) = rho(i,j,k) * vf(i,j,k,n);
                    });
                }
            }

            int face_comp = 0;
            int ncomp = AMREX_SPACEDIM;
            bool is_velocity = true;
            bool allow_inflow_on_outflow = false;
            Array4<int const> const& velbc_arr = velBC_MF ? (*velBC_MF).const_array(mfi)
                                                          : Array4<int const>{};
            HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi,
                                                     (m_advect_momentum) ? rhovel[lev].array(mfi) : vel[lev]->const_array(mfi),
                                                     vel_nph.const_array(mfi),
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
                                                     m_advection_type, PPM::default_limiter,
                                                     allow_inflow_on_outflow, velbc_arr);

            // ************************************************************************
            // Density
            // ************************************************************************
            if (!m_constant_density)
            {
                face_comp = AMREX_SPACEDIM;
                ncomp = 1;
                is_velocity = false;
                allow_inflow_on_outflow = false;
                Array4<int const> const& densbc_arr = densBC_MF ? (*densBC_MF).const_array(mfi)
                                                                : Array4<int const>{};
                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi,
                                                         density[lev]->const_array(mfi),
                                                         rho_nph.const_array(mfi),
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
                                                         m_advection_type, PPM::default_limiter,
                                                         allow_inflow_on_outflow, densbc_arr);
            }

            // ************************************************************************
            // Tracer
            // ************************************************************************
            // Make a FAB holding (rho * tracer) that is the same size as the original tracer FAB
            if (m_advect_tracer && (m_ntrac>0)) {

                // Note we must actually grow the tilebox, not use growntilebox, because
                // we want to use this immediately below and we need all the "ghost cells" of
                // the tiled region
                Box const& bxg = amrex::grow(bx,tracer[lev]->nGrow());

                Array4<Real const> tra     =  tracer[lev]->const_array(mfi);
                Array4<Real const> rho     = density[lev]->const_array(mfi);
                Array4<Real      > trac_tmp;

                auto const* iconserv = get_tracer_iconserv_device_ptr();
                if ( any_conserv_trac ) {
                    trac_tmp = rhotrac[lev].array(mfi);

                    ParallelFor(bxg, m_ntrac,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if ( iconserv[n] ){
                            trac_tmp(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
                        } else {
                            trac_tmp(i,j,k,n) = tra(i,j,k,n);
                        }
                    });
                }

                if (m_constant_density)
                   face_comp = AMREX_SPACEDIM;
                else
                   face_comp = AMREX_SPACEDIM+1;
                ncomp = m_ntrac;
                is_velocity = false;
                allow_inflow_on_outflow = false;
                Array4<int const> const& tracbc_arr = tracBC_MF ? (*tracBC_MF).const_array(mfi)
                                                                : Array4<int const>{};
                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi,
                                          any_conserv_trac ? trac_tmp : tracer[lev]->const_array(mfi),
                                          trac_nph.const_array(mfi),
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
                                          m_advection_type, PPM::default_limiter,
                                          allow_inflow_on_outflow, tracbc_arr);
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
        MultiFab dvdt_tmp(vel[lev]->boxArray(),dmap[lev],AMREX_SPACEDIM,3,MFInfo(),Factory(lev));
        MultiFab drdt_tmp(vel[lev]->boxArray(),dmap[lev],1             ,3,MFInfo(),Factory(lev));
        MultiFab dtdt_tmp(vel[lev]->boxArray(),dmap[lev],m_ntrac       ,3,MFInfo(),Factory(lev));

        // Must initialize to zero because not all values may be set, e.g. outside the domain.
        dvdt_tmp.setVal(0.);
        drdt_tmp.setVal(0.);
        dtdt_tmp.setVal(0.);

        const EBFArrayBoxFactory* ebfact = &EBFactory(lev);
        auto const& vfrac = ebfact->getVolFrac();
#endif

        Real mult = -1.0;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*conv_u[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();

            int flux_comp = 0;
            int  num_comp = AMREX_SPACEDIM;
#ifdef AMREX_USE_EB
            FArrayBox rhovel_fab;
            if (m_advect_momentum && m_eb_flow.enabled)
            {
                // FIXME - not sure if rhovel needs same num grow cells as vel or if
                // could make do with tmp_ng
                Box const& bxg = amrex::grow(bx,vel[lev]->nGrow());
                rhovel_fab.resize(bxg, AMREX_SPACEDIM, The_Async_Arena());

                Array4<Real const> U       = get_velocity_eb()[lev]->const_array(mfi);
                Array4<Real const> rho     =  get_density_eb()[lev]->const_array(mfi);
                Array4<Real      > rho_vel =  rhovel_fab.array();

                ParallelFor(bxg, AMREX_SPACEDIM,
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
                                                 m_eb_flow.enabled ?
                                                    get_velocity_eb()[lev]->const_array(mfi) : Array4<Real const>{},
                                                 m_eb_flow.enabled ?
                                                 ( m_advect_momentum  ? rhovel_fab.const_array() : get_velocity_eb()[lev]->const_array(mfi))
                                                 : Array4<Real const>{},
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
//Was this OMP intentionally left off?
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(*conv_r[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
            EBCellFlagFab const& flagfab = ebfact->getMultiEBCellFlagFab()[mfi];
            if (flagfab.getType(bx) != FabType::covered)
                HydroUtils::EB_ComputeDivergence(bx, drdt_tmp.array(mfi),
                                                 AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                              flux_y[lev].const_array(mfi,flux_comp),
                                                              flux_z[lev].const_array(mfi,flux_comp)),
                                                 vfrac.const_array(mfi), 1, geom[lev], mult,
                                                 fluxes_are_area_weighted,
                                                 m_eb_flow.enabled ?
                                                    get_velocity_eb()[lev]->const_array(mfi) : Array4<Real const>{},
                                                 m_eb_flow.enabled ?
                                                    get_density_eb()[lev]->const_array(mfi) : Array4<Real const>{},
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

        if (m_advect_tracer && m_ntrac > 0)
        {
          int flux_comp = (m_constant_density) ? AMREX_SPACEDIM : AMREX_SPACEDIM+1;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
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
                                                 m_eb_flow.enabled ?
                                                    get_velocity_eb()[lev]->const_array(mfi) : Array4<Real const>{},
                                                 m_eb_flow.enabled ?
                                                    get_tracer_eb()[lev]->const_array(mfi) : Array4<Real const>{},
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

            if ( any_convective_trac )
            {
                // If convective, we define u dot grad trac = div (u trac) - trac div(u)
                HydroUtils::ComputeConvectiveTerm(bx, m_ntrac, mfi,
                                                  tracer[lev]->array(mfi,0),
                                                  AMREX_D_DECL(face_x[lev].array(mfi,flux_comp),
                                                               face_y[lev].array(mfi,flux_comp),
                                                               face_z[lev].array(mfi,flux_comp)),
                                                  divu[lev].array(mfi),
                                                  update_arr,
                                                  get_tracer_iconserv_device_ptr(),
#ifdef AMREX_USE_EB
                                                  *ebfact,
#endif
                                                  m_advection_type);
            }
          } // mfi
        } // advect tracer

#ifdef AMREX_USE_EB
        // We only filled these on the valid cells so we fill same-level interior ghost cells here.
        // (We don't need values outside the domain or at a coarser level so we can call just FillBoundary)
        dvdt_tmp.FillBoundary(geom[lev].periodicity());
        drdt_tmp.FillBoundary(geom[lev].periodicity());
        dtdt_tmp.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*density[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();

            // velocity
            auto const& bc_vel = get_velocity_bcrec_device_ptr();
            redistribute_term(mfi, *conv_u[lev], dvdt_tmp,
                              (m_advect_momentum) ? rhovel[lev] : *vel[lev],
                              bc_vel, lev);

            // density
            if (!m_constant_density) {
                auto const& bc_den = get_density_bcrec_device_ptr();
                redistribute_term(mfi, *conv_r[lev], drdt_tmp,
                                  *density[lev], bc_den, lev);
            } else {
                auto const& drdt = conv_r[lev]->array(mfi);
                ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    drdt(i,j,k) = 0.;
                });
            }

            if (m_advect_tracer) {
                auto const& bc_tra = get_tracer_bcrec_device_ptr();
                redistribute_term(mfi, *conv_t[lev], dtdt_tmp,
                                  any_conserv_trac ? rhotrac[lev] : *tracer[lev],
                                  bc_tra, lev);
            }
        } // mfi
#endif
    } // lev
}
