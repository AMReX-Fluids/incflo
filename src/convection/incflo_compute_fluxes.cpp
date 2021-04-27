#include <Convection.H>
#include <Godunov.H>
#include <MOL.H>

#ifdef AMREX_USE_EB
#include <EBGodunov.H>
#endif

using namespace amrex;

void
convection::compute_fluxes (Box const& bx, MFIter const& mfi,
                            Array4<Real const> const& vel,
                            Array4<Real const> const& rho,
                            Array4<Real const> const& rhotra,
                            Array4<Real const> const& divu,
                            AMREX_D_DECL(Array4<Real const> const& umac,
                                         Array4<Real const> const& vmac,
                                         Array4<Real const> const& wmac),
                            AMREX_D_DECL(Array4<Real> const& fx,
                                         Array4<Real> const& fy,
                                         Array4<Real> const& fz),
                            Array4<Real const> const& fvel,
                            Array4<Real const> const& ftra,
                            Vector<BCRec> const& l_bcrec_velocity,
                                    BCRec const* l_bcrec_velocity_d,
                                      int const* l_iconserv_velocity_d,
                            Vector<BCRec> const& l_bcrec_density,
                                    BCRec const* l_bcrec_density_d,
                                      int const* l_iconserv_density_d,
                            Vector<BCRec> const& l_bcrec_tracer,
                                    BCRec const* l_bcrec_tracer_d,
                                      int const* l_iconserv_tracer_d,
                            std::string l_advection_type, bool l_constant_density,
                            bool l_advect_tracer, int l_ntrac,
                            bool l_godunov_ppm, bool l_godunov_use_forces_in_trans,
#ifdef AMREX_USE_EB
                            EBFArrayBoxFactory const* ebfact,
#endif
                            Geometry& geom, Real l_dt)
{
#ifdef AMREX_USE_EB
    EBCellFlagFab const& flagfab = ebfact->getMultiEBCellFlagFab()[mfi];
    Array4<EBCellFlag const> const& flag = flagfab.const_array();

    bool regular = (flagfab.getType(amrex::grow(bx,2)) == FabType::regular);

    Array4<Real const> AMREX_D_DECL(fcx, fcy, fcz), ccc, vfrac, AMREX_D_DECL(apx, apy, apz);
    if (!regular) {
        AMREX_D_TERM(fcx = ebfact->getFaceCent()[0]->const_array(mfi);,
                     fcy = ebfact->getFaceCent()[1]->const_array(mfi);,
                     fcz = ebfact->getFaceCent()[2]->const_array(mfi););
        ccc   = ebfact->getCentroid().const_array(mfi);
        AMREX_D_TERM(apx = ebfact->getAreaFrac()[0]->const_array(mfi);,
                     apy = ebfact->getAreaFrac()[1]->const_array(mfi);,
                     apz = ebfact->getAreaFrac()[2]->const_array(mfi););
        vfrac = ebfact->getVolFrac().const_array(mfi);
    }
#endif
  
    int nmaxcomp = AMREX_SPACEDIM;
    if (l_advect_tracer) nmaxcomp = std::max(nmaxcomp,l_ntrac);

        int n_tmp_fac;
#if (AMREX_SPACEDIM == 3)
        n_tmp_fac = 14;
#else
        n_tmp_fac = 10;
#endif

        // n_tmp_grow is used to create tmpfab which is passed in to compute_godunov_fluxes
        // (both regular and EB) as pointer "p" and is used to hold Imx/Ipx etc ...
        int n_tmp_grow;
#ifdef AMREX_USE_EB
        n_tmp_grow = 4;
#else
        n_tmp_grow = 4;
#endif

        FArrayBox tmpfab(amrex::grow(bx,n_tmp_grow), nmaxcomp*n_tmp_fac+1);
        Elixir eli = tmpfab.elixir();

#ifdef AMREX_USE_EB
        if (!regular)
        {
            int flux_comp = 0;
            if (l_advection_type == "Godunov")
                ebgodunov::compute_godunov_fluxes(bx, flux_comp, AMREX_SPACEDIM,
                                                  AMREX_D_DECL(fx, fy, fz), vel, 
                                                  AMREX_D_DECL(umac, vmac, wmac), 
                                                  fvel, divu, l_dt, 
                                                  l_bcrec_velocity,
                                                  l_bcrec_velocity_d,
                                                  l_iconserv_velocity_d,
                                                  tmpfab.dataPtr(), flag, 
                                                  AMREX_D_DECL(apx, apy, apz), vfrac,
                                                  AMREX_D_DECL(fcx, fcy, fcz), ccc, 
                                                  geom, true); // is_velocity
            else
                MOL::compute_convective_fluxes_eb(bx, flux_comp, AMREX_SPACEDIM,
                                                  AMREX_D_DECL(fx, fy, fz), vel, 
                                                  AMREX_D_DECL(umac, vmac, wmac),
                                                  l_bcrec_velocity.data(),
                                                  l_bcrec_velocity_d,
                                                  l_iconserv_velocity_d,
                                                  flag, AMREX_D_DECL(fcx, fcy, fcz), ccc, vfrac, geom);
            flux_comp += AMREX_SPACEDIM;

            if (!l_constant_density) {
                if (l_advection_type == "Godunov")
                    ebgodunov::compute_godunov_fluxes(bx, flux_comp, 1,
                                                      AMREX_D_DECL(fx, fy, fz), rho,
                                                      AMREX_D_DECL(umac, vmac, wmac), 
                                                      {}, divu, l_dt, 
                                                      l_bcrec_density,
                                                      l_bcrec_density_d,
                                                      l_iconserv_density_d,
                                                      tmpfab.dataPtr(), flag,
                                                      AMREX_D_DECL(apx, apy, apz), vfrac,
                                                      AMREX_D_DECL(fcx, fcy, fcz), ccc, 
                                                      geom);
                else
                    MOL::compute_convective_fluxes_eb(bx, flux_comp, 1,
                                                      AMREX_D_DECL(fx, fy, fz), rho, 
                                                      AMREX_D_DECL(umac, vmac, wmac),
                                                      l_bcrec_density.data(),
                                                      l_bcrec_density_d,
                                                      l_iconserv_density_d,
                                                      flag, AMREX_D_DECL(fcx, fcy, fcz), ccc, vfrac, geom);
                flux_comp += 1;
            }
            if (l_advect_tracer) {
                if (l_advection_type == "Godunov")
                    ebgodunov::compute_godunov_fluxes(bx, flux_comp, l_ntrac,
                                                      AMREX_D_DECL(fx, fy, fz), rhotra,
                                                      AMREX_D_DECL(umac, vmac, wmac), 
                                                      ftra, divu, l_dt, 
                                                      l_bcrec_tracer,
                                                      l_bcrec_tracer_d,
                                                      l_iconserv_tracer_d,
                                                      tmpfab.dataPtr(), flag,
                                                      AMREX_D_DECL(apx, apy, apz), vfrac,
                                                      AMREX_D_DECL(fcx, fcy, fcz), ccc, 
                                                      geom);
                else
                    MOL::compute_convective_fluxes_eb(bx, flux_comp, l_ntrac,
                                                      AMREX_D_DECL(fx, fy, fz), rhotra, 
                                                      AMREX_D_DECL(umac, vmac, wmac),
                                                      l_bcrec_tracer.data(),
                                                      l_bcrec_tracer_d,
                                                      l_iconserv_tracer_d,
                                                      flag, AMREX_D_DECL(fcx, fcy, fcz), ccc, vfrac, geom);
            }
            Gpu::streamSynchronize();
        }
        else
#endif
        {
            int flux_comp = 0;
                if (l_advection_type == "Godunov")
                    godunov::compute_godunov_fluxes(bx, flux_comp, AMREX_SPACEDIM,
                                                    AMREX_D_DECL(fx, fy, fz), vel, 
                                                    AMREX_D_DECL(umac, vmac, wmac), 
                                                    fvel, divu, l_dt, 
                                                    l_bcrec_velocity_d,
                                                    l_iconserv_velocity_d,
                                                    tmpfab.dataPtr(),l_godunov_ppm, 
                                                    l_godunov_use_forces_in_trans, 
                                                    geom, true);
                else
                    MOL::compute_convective_fluxes(bx, flux_comp, AMREX_SPACEDIM, AMREX_D_DECL(fx, fy, fz), vel,
                                                   AMREX_D_DECL(umac, vmac, wmac),
                                                   l_bcrec_velocity.data(),
                                                   l_bcrec_velocity_d,
                                                   l_iconserv_velocity_d, geom);

            flux_comp += AMREX_SPACEDIM;

            if (!l_constant_density) {
                if (l_advection_type == "Godunov")
                    godunov::compute_godunov_fluxes(bx, flux_comp, 1,
                                                    AMREX_D_DECL(fx, fy, fz), rho, 
                                                    AMREX_D_DECL(umac, vmac, wmac),
                                                    {}, divu, l_dt, 
                                                    l_bcrec_density_d,
                                                    l_iconserv_density_d,
                                                    tmpfab.dataPtr(),l_godunov_ppm,
                                                    l_godunov_use_forces_in_trans,
                                                    geom);
                else
                    MOL::compute_convective_fluxes(bx, flux_comp, 1, AMREX_D_DECL(fx, fy, fz), rho,
                                                   AMREX_D_DECL(umac, vmac, wmac),
                                                   l_bcrec_density.data(),
                                                   l_bcrec_density_d, 
                                                   l_iconserv_density_d, geom);
               flux_comp += 1;
            }
            if (l_advect_tracer) {
                if (l_advection_type == "Godunov")
                    godunov::compute_godunov_fluxes(bx, flux_comp, l_ntrac,
                                                    AMREX_D_DECL(fx, fy, fz), rhotra, 
                                                    AMREX_D_DECL(umac, vmac, wmac),
                                                    ftra, divu, l_dt, 
                                                    l_bcrec_tracer_d,
                                                    l_iconserv_tracer_d,
                                                    tmpfab.dataPtr(),l_godunov_ppm,
                                                    l_godunov_use_forces_in_trans,
                                                    geom);
                else
                   MOL::compute_convective_fluxes(bx, flux_comp, l_ntrac, AMREX_D_DECL(fx, fy, fz), rhotra,
                                                  AMREX_D_DECL(umac, vmac, wmac),
                                                  l_bcrec_tracer.data(),
                                                  l_bcrec_tracer_d, 
                                                  l_iconserv_tracer_d, geom);
            }
            Gpu::streamSynchronize();
        }
}
