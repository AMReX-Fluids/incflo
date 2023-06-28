/**
 * \file hydro_redistribution.H
 * \addtogroup Redistribution
 * @{
 *
 */

#ifndef REDISTRIBUTION_H_
#define REDISTRIBUTION_H_

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiCutFab.H>

/**
 * Placeholder description of Redistribution namespace.
 *
 */

namespace Redistribution {

    namespace defaults {
        constexpr int srd_max_order = 2;
        constexpr amrex::Real target_vol_fraction = 0.5;
        constexpr amrex::Real small_norm_diff = 1.e-8;
    }

     void FillNewlyUncovered ( amrex::MultiFab& mf,
                               amrex::EBFArrayBoxFactory const& ebfact_old,
                               amrex::EBFArrayBoxFactory const& ebfact_new,
                               amrex::MultiFab const& vel_eb,
                               amrex::Geometry& geom,
                               amrex::Real target_volfrac);

    // FIXME - probably would benefit from passing ebfact instead of all the pieces...
// Everything should eventually come down to calling this function
    void Apply ( amrex::Box const& bx, int ncomp,
                 amrex::Array4<amrex::Real>       const& out,
                 amrex::Array4<amrex::Real>       const& dUdt_in,
                 amrex::Array4<amrex::Real const> const& U_in,
                 amrex::Array4<amrex::Real> const& scratch,
                 amrex::Array4<amrex::EBCellFlag const> const& flag_old,
                 amrex::Array4<amrex::EBCellFlag const> const& flag_new,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx_old,
                              amrex::Array4<amrex::Real const> const& apy_old,
                              amrex::Array4<amrex::Real const> const& apz_old),
                 amrex::Array4<amrex::Real const> const& vfrac_old,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx_new,
                              amrex::Array4<amrex::Real const> const& apy_new,
                              amrex::Array4<amrex::Real const> const& apz_new),
                 amrex::Array4<amrex::Real const> const& vfrac_new,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                              amrex::Array4<amrex::Real const> const& fcy,
                              amrex::Array4<amrex::Real const> const& fcz),
                 amrex::Array4<amrex::Real const> const& ccent,
                 amrex::BCRec  const* d_bcrec_ptr,
                 amrex::Geometry const& geom,
                 amrex::Real dt, std::string redistribution_type,
                 amrex::Array4<amrex::Real const> const& vel_eb_old,
                 amrex::Array4<amrex::Real const> const& bnorm_old,
                 amrex::Array4<amrex::Real const> const& barea_old,
                 amrex::Array4<amrex::Real const> const& vel_eb_new,
                 amrex::Array4<amrex::Real const> const& bnorm_new,
                 amrex::Array4<amrex::Real const> const& barea_new,
                 const int max_order,
                 amrex::Real target_volfrac,
                 amrex::Array4<amrex::Real const> const& update_scale);

    // for backwards compatibility
    void Apply ( amrex::Box const& bx, int ncomp,
                 amrex::Array4<amrex::Real>       const& out,
                 amrex::Array4<amrex::Real>       const& dUdt_in,
                 amrex::Array4<amrex::Real const> const& U_in,
                 amrex::Array4<amrex::Real> const& scratch,
                 amrex::Array4<amrex::EBCellFlag const> const& flag,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                              amrex::Array4<amrex::Real const> const& apy,
                              amrex::Array4<amrex::Real const> const& apz),
                 amrex::Array4<amrex::Real const> const& vfrac,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                              amrex::Array4<amrex::Real const> const& fcy,
                              amrex::Array4<amrex::Real const> const& fcz),
                 amrex::Array4<amrex::Real const> const& ccent,
                 amrex::BCRec  const* d_bcrec_ptr,
                 amrex::Geometry const& geom,
                 amrex::Real dt, std::string redistribution_type,
                 const int max_order = defaults::srd_max_order,
                 amrex::Real target_volfrac = defaults::target_vol_fraction,
                 amrex::Array4<amrex::Real const> const& update_scale={});

    // FIXME? This initial one is always stationary right? then we can remove a vfrac...
    void ApplyToInitialData ( amrex::Box const& bx, int ncomp,
                              amrex::Array4<amrex::Real                  > const& U_out,
                              amrex::Array4<amrex::Real                  > const& U_in,
                              amrex::Array4<amrex::EBCellFlag const> const& flag,
                              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                           amrex::Array4<amrex::Real const> const& apy,
                                           amrex::Array4<amrex::Real const> const& apz),
                              amrex::Array4<amrex::Real const> const& vfrac_old,
                              amrex::Array4<amrex::Real const> const& vfrac_new,
                              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                           amrex::Array4<amrex::Real const> const& fcy,
                                           amrex::Array4<amrex::Real const> const& fcz),
                              amrex::Array4<amrex::Real const> const& ccent,
                              amrex::BCRec  const* d_bcrec_ptr,
                              amrex::Geometry& geom, std::string redistribution_type,
                              const int srd_max_order = defaults::srd_max_order,
                              amrex::Real target_volfrac = defaults::target_vol_fraction);

    void StateRedistribute ( amrex::Box const& bx, int ncomp,
                             amrex::Array4<amrex::Real> const& dUdt_out,
                             amrex::Array4<amrex::Real> const& dUdt_in,
                             amrex::Array4<amrex::EBCellFlag const> const& flag_new,
                             amrex::Array4<amrex::Real const> const& vfrac_old,
                             amrex::Array4<amrex::Real const> const& vfrac_new,
                             AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                          amrex::Array4<amrex::Real const> const& fcy,
                                          amrex::Array4<amrex::Real const> const& fcz),
                             amrex::Array4<amrex::Real const> const& ccent,
                             amrex::BCRec  const* d_bcrec_ptr,
                             amrex::Array4<int const> const& itracker,
                             amrex::Array4<amrex::Real const> const& nrs,
                             amrex::Array4<amrex::Real const> const& alpha,
                             amrex::Array4<amrex::Real const> const& nbhd_vol,
                             amrex::Array4<amrex::Real const> const& cent_hat,
                             amrex::Geometry const& geom,
                             const int max_order = defaults::srd_max_order);

    void MakeITracker ( amrex::Box const& bx,
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx_old,
                                     amrex::Array4<amrex::Real const> const& apy_old,
                                     amrex::Array4<amrex::Real const> const& apz_old),
                        amrex::Array4<amrex::Real const> const& vfrac_old,
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx_new,
                                     amrex::Array4<amrex::Real const> const& apy_new,
                                     amrex::Array4<amrex::Real const> const& apz_new),
                        amrex::Array4<amrex::Real const> const& vfrac_new,
                        amrex::Array4<int> const& itracker,
                        amrex::Geometry const& geom,
                        amrex::Real target_volfrac,
                        amrex::Array4<amrex::Real const> const& vel_eb = {});

    void MakeStateRedistUtils ( amrex::Box const& bx,
                                amrex::Array4<amrex::Real const> const& vfrac_old,
                                amrex::Array4<amrex::Real const> const& vfrac_new,
                                amrex::Array4<amrex::Real const> const& ccent,
                                amrex::Array4<        int const> const& itracker,
                                amrex::Array4<amrex::Real> const& nrs,
                                amrex::Array4<amrex::Real> const& nbhd_vol,
                                amrex::Array4<amrex::Real> const& alpha,
                                amrex::Array4<amrex::Real> const& cent_hat,
                                amrex::Geometry const& geom,
                                amrex::Real target_volfrac);

} // namespace redistribution

#endif
/** @} */
