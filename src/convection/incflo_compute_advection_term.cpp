#include <incflo.H>
#include <hydro_godunov.H>
#include <hydro_mol.H>
#include <hydro_utils.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#include <hydro_ebgodunov.H>
#include <hydro_ebmol.H>
#endif

using namespace amrex;

void incflo::init_advection ()
{
#ifdef AMREX_USE_EB
    // We default to conservative if using EB
    m_iconserv_velocity.resize(AMREX_SPACEDIM, 1);
    m_iconserv_velocity_d.resize(AMREX_SPACEDIM, 1);
#else
    // We default to non-conservative if not using EB
    m_iconserv_velocity.resize(AMREX_SPACEDIM, 0);
    m_iconserv_velocity_d.resize(AMREX_SPACEDIM, 0);
#endif

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
    int ngmac = nghost_mac();

    bool fluxes_are_area_weighted = false;

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

    Vector<Array<MultiFab*,AMREX_SPACEDIM> > fluxes(finest_level+1);
    Vector<Array<MultiFab*,AMREX_SPACEDIM> >  faces(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(
           face_x[lev].define(u_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),Factory(lev));,
           face_y[lev].define(v_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),Factory(lev));,
           face_z[lev].define(w_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),Factory(lev)););
        AMREX_D_TERM(
           flux_x[lev].define(u_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),Factory(lev));,
           flux_y[lev].define(v_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),Factory(lev));,
           flux_z[lev].define(w_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),Factory(lev)););

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

    int ngrow = 4;

#ifdef AMREX_USE_EB
    if (m_redistribution_type=="StateRedist")
        ++ngrow;
#endif

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (ngmac > 0) {
            AMREX_D_TERM(u_mac[lev]->FillBoundary(geom[lev].periodicity());,
                         v_mac[lev]->FillBoundary(geom[lev].periodicity());,
                         w_mac[lev]->FillBoundary(geom[lev].periodicity()););
        }

        MultiFab divu(vel[lev]->boxArray(),vel[lev]->DistributionMap(),1,4);
        divu.setVal(0.);
        Array<MultiFab const*, AMREX_SPACEDIM> u;
        AMREX_D_TERM(u[0] = u_mac[lev];,
                     u[1] = v_mac[lev];,
                     u[2] = w_mac[lev];);
                     
#ifdef AMREX_USE_EB
        Array4<Real const> AMREX_D_DECL(fcx, fcy, fcz), ccc, vfrac, AMREX_D_DECL(apx, apy, apz);

        const auto& ebfact = EBFactory(lev);

        if (!ebfact.isAllRegular())
            amrex::EB_computeDivergence(divu,u,geom[lev],true);
        else
#endif
        amrex::computeDivergence(divu,u,geom[lev]);

        divu.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*density[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();

            Array4<Real const> const& divu_arr = divu.const_array(mfi);

#ifdef AMREX_USE_EB
            EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
            Array4<EBCellFlag const> const& flag = flagfab.const_array();

            bool regular = (flagfab.getType(amrex::grow(bx,2)) == FabType::regular);

            if (!regular)
            {
                vfrac = ebfact.getVolFrac().const_array(mfi);
                ccc   = ebfact.getCentroid().const_array(mfi);

                AMREX_D_TERM( apx = ebfact.getAreaFrac()[0]->const_array(mfi);,
                              apy = ebfact.getAreaFrac()[1]->const_array(mfi);,
                              apz = ebfact.getAreaFrac()[2]->const_array(mfi););

                AMREX_D_TERM( fcx = ebfact.getFaceCent()[0]->const_array(mfi);,
                              fcy = ebfact.getFaceCent()[1]->const_array(mfi);,
                              fcz = ebfact.getFaceCent()[2]->const_array(mfi););
            }
#endif
            // ************************************************************************
            // Velocity
            // ************************************************************************
            int face_comp = 0;  
            int ncomp = AMREX_SPACEDIM;
            bool is_velocity = true;
            HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi, 
                                                  vel[lev]->const_array(mfi),
                                     AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                  flux_y[lev].array(mfi,face_comp),
                                                  flux_z[lev].array(mfi,face_comp)),
                                     AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                  face_y[lev].array(mfi,face_comp),
                                                  face_z[lev].array(mfi,face_comp)),
                                     AMREX_D_DECL(u_mac[lev]->const_array(mfi),
                                                  v_mac[lev]->const_array(mfi),
                                                  w_mac[lev]->const_array(mfi)),
                                     divu_arr,
                                     (!vel_forces.empty()) ? vel_forces[lev]->const_array(mfi) : Array4<Real const>{},
                                     geom[lev], m_dt, 
                                     get_velocity_bcrec(), 
                                     get_velocity_bcrec_device_ptr(),
                                     get_velocity_iconserv_device_ptr(),
#ifdef AMREX_USE_EB
                                     ebfact,
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
#endif
                                          m_godunov_ppm, m_godunov_use_forces_in_trans,
                                          is_velocity, fluxes_are_area_weighted,
                                          m_advection_type);
            }

            // ************************************************************************
            // (Rho*Tracer)
            // ************************************************************************
            // Make a FAB holding (rho * tracer) that is the same size as the original tracer FAB
            FArrayBox rhotracfab;
            if (m_advect_tracer && (m_ntrac>0)) {

                Box rhotrac_box = Box((*tracer[lev])[mfi].box());
                Elixir eli_rt;
                Array4<Real> rhotrac;
                Array4<Real const> tra =  tracer[lev]->const_array(mfi);
                Array4<Real const> rho = density[lev]->const_array(mfi);
                rhotracfab.resize(rhotrac_box, m_ntrac);
                eli_rt  = rhotracfab.elixir();
                rhotrac = rhotracfab.array();
                amrex::ParallelFor(rhotrac_box, m_ntrac,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhotrac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
                });

                if (m_constant_density)
                   face_comp = AMREX_SPACEDIM;
                else
                   face_comp = AMREX_SPACEDIM+1;
                ncomp = m_ntrac;
                is_velocity = false;

                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi, 
                                                       rhotrac,
                                          AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                       flux_y[lev].array(mfi,face_comp),
                                                       flux_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                       face_y[lev].array(mfi,face_comp),
                                                       face_z[lev].array(mfi,face_comp)),
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

    int flux_comp;

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

#ifdef AMREX_USE_EB
//      if (!ebfact->isAllRegular())
//          amrex::EB_computeDivergence(divu,u,geom[lev],true);
//      else
#endif
//      amrex::computeDivergence(divu,u,geom[lev]);

        Real mult = -1.0;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*conv_u[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();

            flux_comp = 0;
#ifdef AMREX_USE_EB
            EBCellFlagFab const& flagfab = ebfact->getMultiEBCellFlagFab()[mfi];
            if (flagfab.getType(bx) != FabType::covered)
                HydroUtils::EB_ComputeDivergence(bx, dvdt_tmp.array(mfi),
                                                 AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                              flux_y[lev].const_array(mfi,flux_comp),
                                                              flux_z[lev].const_array(mfi,flux_comp)),
                                                 vfrac.const_array(mfi), AMREX_SPACEDIM, geom[lev],
                                                 mult,
                                                 fluxes_are_area_weighted);
#else
            HydroUtils::ComputeDivergence(bx, conv_u[lev]->array(mfi),
                                          AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                       flux_y[lev].const_array(mfi,flux_comp),
                                                       flux_z[lev].const_array(mfi,flux_comp)),
                                          AMREX_D_DECL(face_x[lev].const_array(mfi,flux_comp),
                                                       face_y[lev].const_array(mfi,flux_comp),
                                                       face_z[lev].const_array(mfi,flux_comp)),
                                          AMREX_D_DECL(u_mac[lev]->const_array(mfi),
                                                       v_mac[lev]->const_array(mfi),
                                                       w_mac[lev]->const_array(mfi)),
                                          AMREX_SPACEDIM, geom[lev],
                                          get_velocity_iconserv_device_ptr(), mult,
                                          fluxes_are_area_weighted);
#endif
        } // mfi

        if (!m_constant_density)
        {
          flux_comp = AMREX_SPACEDIM;
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
                                                 fluxes_are_area_weighted);
#else
            HydroUtils::ComputeDivergence(bx, conv_r[lev]->array(mfi),
                                          AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                       flux_y[lev].const_array(mfi,flux_comp),
                                                       flux_z[lev].const_array(mfi,flux_comp)),
                                          AMREX_D_DECL(face_x[lev].const_array(mfi,flux_comp),
                                                       face_y[lev].const_array(mfi,flux_comp),
                                                       face_z[lev].const_array(mfi,flux_comp)),
                                          AMREX_D_DECL(u_mac[lev]->const_array(mfi),
                                                       v_mac[lev]->const_array(mfi),
                                                       w_mac[lev]->const_array(mfi)),
                                          1, geom[lev],
                                          get_density_iconserv_device_ptr(), mult,
                                          fluxes_are_area_weighted);
#endif
          } // mfi
        } // not constant density

        if (m_advect_tracer && m_ntrac > 0)
        {
          if (m_constant_density)
              flux_comp = AMREX_SPACEDIM;
          else
              flux_comp = AMREX_SPACEDIM+1;
          for (MFIter mfi(*conv_t[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
            EBCellFlagFab const& flagfab = ebfact->getMultiEBCellFlagFab()[mfi];
            if (flagfab.getType(bx) != FabType::covered)
                HydroUtils::EB_ComputeDivergence(bx, dtdt_tmp.array(mfi),
                                                 AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                              flux_y[lev].const_array(mfi,flux_comp),
                                                              flux_z[lev].const_array(mfi,flux_comp)),
                                                 vfrac.const_array(mfi), m_ntrac, geom[lev], mult,
                                                 fluxes_are_area_weighted);
#else
            HydroUtils::ComputeDivergence(bx, conv_t[lev]->array(mfi),
                                          AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                       flux_y[lev].const_array(mfi,flux_comp),
                                                       flux_z[lev].const_array(mfi,flux_comp)),
                                          AMREX_D_DECL(face_x[lev].const_array(mfi,flux_comp),
                                                       face_y[lev].const_array(mfi,flux_comp),
                                                       face_z[lev].const_array(mfi,flux_comp)),
                                          AMREX_D_DECL(u_mac[lev]->const_array(mfi),
                                                       v_mac[lev]->const_array(mfi),
                                                       w_mac[lev]->const_array(mfi)),
                                          m_ntrac, geom[lev],
                                          get_tracer_iconserv_device_ptr(), mult,
                                          fluxes_are_area_weighted);
#endif
          } // mfi
        } // advect tracer

#ifdef AMREX_USE_EB
        // We only filled these on the valid cells so we fill same-level interior ghost cells here. 
        // (We don't need values outside the domain or at a coarser level so we can call just FillBoundary)
        dvdt_tmp.FillBoundary(geom[lev].periodicity());
        drdt_tmp.FillBoundary(geom[lev].periodicity());
        dtdt_tmp.FillBoundary(geom[lev].periodicity());

        for (MFIter mfi(*density[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Make a FAB holding (rho * tracer) that is the same size as the original tracer FAB
            FArrayBox rhotracfab;
            if (m_advect_tracer && (m_ntrac > 0)) {
                Box rhotrac_box = Box((*tracer[lev])[mfi].box());
                Elixir eli_rt;
                Array4<Real> rhotrac;
                Array4<Real const> tra =  tracer[lev]->const_array(mfi);
                Array4<Real const> rho = density[lev]->const_array(mfi);
                rhotracfab.resize(rhotrac_box, m_ntrac);
                eli_rt  = rhotracfab.elixir();
                rhotrac = rhotracfab.array();
                amrex::ParallelFor(rhotrac_box, m_ntrac,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhotrac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
                });
            }

            Box const& bx = mfi.tilebox();
            redistribute_convective_term (bx, mfi,
                                          vel[lev]->const_array(mfi),
                                          density[lev]->const_array(mfi),
                                          (m_advect_tracer && (m_ntrac>0)) ? rhotracfab.const_array() : Array4<Real const>{},
                                          dvdt_tmp.array(mfi),
                                          drdt_tmp.array(mfi),
                                          (m_ntrac>0) ? dtdt_tmp.array(mfi) : Array4<Real>{},
                                          conv_u[lev]->array(mfi),
                                          conv_r[lev]->array(mfi),
                                          (m_ntrac>0) ? conv_t[lev]->array(mfi) : Array4<Real>{},
                                          m_redistribution_type, m_constant_density, m_advect_tracer, m_ntrac,
                                          ebfact, geom[lev], m_dt);
       } // mfi
#endif
    } // lev
}
