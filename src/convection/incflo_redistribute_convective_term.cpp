#include <AMReX_Config.H>

#ifdef AMREX_USE_EB

#include <hydro_redistribution.H>
//#include <AMReX_MultiFabUtil.H>
//#include <AMReX_MultiCutFab.H>
#include <incflo.H>

using namespace amrex;

void
incflo::redistribute_convective_term ( Box const& bx, MFIter const& mfi,
                                       Array4<Real const > const& vel, // velocity
                                       Array4<Real const > const& rho, // density
                                       Array4<Real const > const& rhotrac, // tracer
                                       Array4<Real> const& dvdt_tmp, // velocity
                                       Array4<Real> const& drdt_tmp, // density
                                       Array4<Real> const& dtdt_tmp, // tracer
                                       Array4<Real> const& dvdt, // velocity
                                       Array4<Real> const& drdt, // density
                                       Array4<Real> const& dtdt, // tracer
                                       std::string l_redistribution_type,
                                       bool l_constant_density,
                                       bool l_advect_tracer, int l_ntrac,
                                       EBFArrayBoxFactory const* ebfact_old,
                                       EBFArrayBoxFactory const* ebfact_new,
                                       Geometry& lev_geom, Real l_dt)
{
    EBCellFlagFab const& flagfab = ebfact_old->getMultiEBCellFlagFab()[mfi];
    Array4<EBCellFlag const> const& flag = flagfab.const_array();

    bool regular = (flagfab.getType(amrex::grow(bx,4)) == FabType::regular);

    Array4<Real const> AMREX_D_DECL(fcx, fcy, fcz);
    Array4<Real const> AMREX_D_DECL(apx_old, apy_old, apz_old);
    Array4<Real const> AMREX_D_DECL(apx_new, apy_new, apz_new);
    Array4<Real const> ccc, vfrac_old, vfrac_new;

    if (!regular)
    {
        AMREX_D_TERM(fcx = ebfact_old->getFaceCent()[0]->const_array(mfi);,
                     fcy = ebfact_old->getFaceCent()[1]->const_array(mfi);,
                     fcz = ebfact_old->getFaceCent()[2]->const_array(mfi););
        ccc   = ebfact_old->getCentroid().const_array(mfi);
        AMREX_D_TERM(apx_old = ebfact_old->getAreaFrac()[0]->const_array(mfi);,
                     apy_old = ebfact_old->getAreaFrac()[1]->const_array(mfi);,
                     apz_old = ebfact_old->getAreaFrac()[2]->const_array(mfi););
        AMREX_D_TERM(apx_new = ebfact_new->getAreaFrac()[0]->const_array(mfi);,
                     apy_new = ebfact_new->getAreaFrac()[1]->const_array(mfi);,
                     apz_new = ebfact_new->getAreaFrac()[2]->const_array(mfi););
        vfrac_old = ebfact_old->getVolFrac().const_array(mfi);
        vfrac_new = ebfact_new->getVolFrac().const_array(mfi);

        Box gbx = bx;

        if (l_redistribution_type == "StateRedist") {
            gbx.grow(3);
        } else if (l_redistribution_type == "FluxRedist") {
            gbx.grow(2);
        }

        int nmaxcomp = AMREX_SPACEDIM;
        if (l_advect_tracer)
            nmaxcomp = std::max(nmaxcomp,l_ntrac);

        FArrayBox scratch_fab(gbx,nmaxcomp);
        Array4<Real> scratch = scratch_fab.array();
        Elixir eli_scratch = scratch_fab.elixir();

        // This is scratch space if calling StateRedistribute
        //  but is used as the weights (here set to 1) if calling
        //  FluxRedistribute
        amrex::ParallelFor(Box(scratch),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                scratch(i,j,k) = 1.;
            });

        // velocity
        auto const& bc_vel = get_velocity_bcrec_device_ptr();
        amrex::Print() << "Redist for velocity " << std::endl;
        // std::string vel_redistribution_type = "NoRedist";

        Redistribution::Apply(bx, AMREX_SPACEDIM, dvdt, dvdt_tmp, vel, scratch, flag,
                              AMREX_D_DECL(apx_old, apy_old, apz_old), vfrac_old,
                              AMREX_D_DECL(apx_new, apy_new, apz_new), vfrac_new,
                              AMREX_D_DECL(fcx, fcy, fcz), ccc,
                              bc_vel, lev_geom, l_dt, l_redistribution_type);

        // density
        if (!l_constant_density) {
            auto const& bc_den = get_density_bcrec_device_ptr();
            amrex::Print() << "Redist for density " << std::endl;
            Redistribution::Apply(bx, 1, drdt, drdt_tmp, rho, scratch, flag,
                                  AMREX_D_DECL(apx_old, apy_old, apz_old), vfrac_old,
                                  AMREX_D_DECL(apx_new, apy_new, apz_new), vfrac_new,
                                  AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                  bc_den, lev_geom, l_dt, l_redistribution_type);
        } else {
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                drdt(i,j,k) = 0.;
            });
        }

        if (l_advect_tracer) {
            auto const& bc_tra = get_tracer_bcrec_device_ptr();
            amrex::Print() << "Redist for tracer " << std::endl;
            Redistribution::Apply(bx, l_ntrac, dtdt, dtdt_tmp, rhotrac, scratch, flag,
                                  AMREX_D_DECL(apx_old, apy_old, apz_old), vfrac_old,
                                  AMREX_D_DECL(apx_new, apy_new, apz_new), vfrac_new,
                                  AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                  bc_tra, lev_geom, l_dt, l_redistribution_type);
        }

    } else {
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int n = 0; n < AMREX_SPACEDIM; n++)
               dvdt(i,j,k,n) = dvdt_tmp(i,j,k,n);

            if (!l_constant_density)
                drdt(i,j,k) = drdt_tmp(i,j,k);

            if (l_advect_tracer)
               for (int n = 0; n < l_ntrac; n++)
                   dtdt(i,j,k,n) = dtdt_tmp(i,j,k,n);
        });
    }
}
#endif
