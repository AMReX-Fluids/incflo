/**
 * \file hydro_redistribution.H
 * \addtogroup Redistribution
 * @{
 *
 */

#ifndef IAMR_REDISTRIBUTION_H_
#define IAMR_REDISTRIBUTION_H_

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiCutFab.H>

/**
 * Placeholder description of Redistribution namespace.
 *
 */

namespace Redistribution {

    namespace defaults {
        constexpr int srd_max_order = 2;
        constexpr amrex::Real target_vol_fraction = 1.0;
    }
    amrex::Array<amrex::Array<int,9>,AMREX_SPACEDIM> getCellMap();
    amrex::Array<int,9> getInvCellMap();

    void Apply ( amrex::Box const& bx, int ncomp,
                 amrex::Array4<amrex::Real>       const& dUdt_out,
                 amrex::Array4<amrex::Real>       const& dUdt_in,
                 amrex::Array4<amrex::Real const> const& U_in,
                 amrex::Array4<amrex::Real> const& scratch,
                 amrex::Array4<amrex::EBCellFlag const> const& flag,
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
                 amrex::Array4<amrex::Real const> const& vel_eb,
                 amrex::Array4<amrex::Real const> const& bnorm,
                 amrex::Array4<amrex::Real const> const& barea,
                 const int max_order,
                 amrex::Real target_volfrac,
                 amrex::Array4<amrex::Real const> const& update_scale);

    // Overloading in this way doens't work. Get multiple defs during linking...
    //
    // void Apply ( amrex::Box const& bx, int ncomp,
    //              amrex::Array4<amrex::Real>       const& dUdt_out,
    //              amrex::Array4<amrex::Real>       const& dUdt_in,
    //              amrex::Array4<amrex::Real const> const& U_in,
    //              amrex::Array4<amrex::Real> const& scratch,
    //              amrex::Array4<amrex::EBCellFlag const> const& flag,
    //              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
    //                           amrex::Array4<amrex::Real const> const& apy,
    //                           amrex::Array4<amrex::Real const> const& apz),
    //              amrex::Array4<amrex::Real const> const& vfrac,
    //              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
    //                           amrex::Array4<amrex::Real const> const& fcy,
    //                           amrex::Array4<amrex::Real const> const& fcz),
    //              amrex::Array4<amrex::Real const> const& ccent,
    //              amrex::BCRec  const* d_bcrec_ptr,
    //              amrex::Geometry const& geom,
    //              amrex::Real dt, std::string redistribution_type,
    //              const int max_order = defaults::srd_max_order,
    //              amrex::Real target_volfrac = defaults::target_vol_fraction,
    //              amrex::Array4<amrex::Real const> const& update_scale={})
    // {
    //  // Apply(bx, ncomp, dUdt_out, dUdt_in, U_in, scratch, flag,
    //  //       AMREX_D_DECL(apx,apy, apz), vfrac,
    //  //       AMREX_D_DECL(amrex::Array4<amrex::Real const> {},
    //  //                 amrex::Array4<amrex::Real const> {},
    //  //                 amrex::Array4<amrex::Real const> {}),
    //  //       amrex::Array4<amrex::Real const> {},
    //  //       AMREX_D_DECL(fcx,fcy, fcz), ccent,
    //  //       d_bcrec_ptr, geom, dt, redistribution_type,
    //  //       amrex::Array4<amrex::Real const> {},
    //  //       amrex::Array4<amrex::Real const> {},
    //  //       amrex::Array4<amrex::Real const> {},
    //  //       max_order, target_volfrac,
    //  //       update_scale);
    // }

    // void Apply ( amrex::Box const& bx, int ncomp,
    //              amrex::Array4<amrex::Real>       const& dUdt_out,
    //              amrex::Array4<amrex::Real>       const& dUdt_in,
    //              amrex::Array4<amrex::Real const> const& U_in,
    //              amrex::Array4<amrex::Real> const& scratch,
    //              amrex::Array4<amrex::EBCellFlag const> const& flag,
    //              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx_old,
    //                           amrex::Array4<amrex::Real const> const& apy_old,
    //                           amrex::Array4<amrex::Real const> const& apz_old),
    //              amrex::Array4<amrex::Real const> const& vfrac_old,
    //              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx_new,
    //                           amrex::Array4<amrex::Real const> const& apy_new,
    //                           amrex::Array4<amrex::Real const> const& apz_new),
    //              amrex::Array4<amrex::Real const> const& vfrac_new,
    //              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
    //                           amrex::Array4<amrex::Real const> const& fcy,
    //                           amrex::Array4<amrex::Real const> const& fcz),
    //              amrex::Array4<amrex::Real const> const& ccent,
    //              amrex::BCRec  const* d_bcrec_ptr,
    //              amrex::Geometry const& geom,
    //              amrex::Real dt, std::string redistribution_type,
    //           amrex::Array4<amrex::Real const> const& vel_eb,
    //           amrex::Array4<amrex::Real const> const& bnorm,
    //           amrex::Array4<amrex::Real const> const& barea)
    // {
    //  // Apply(bx, ncomp, dUdt_out, dUdt_in, U_in, scratch, flag,
    //  //       AMREX_D_DECL(apx_old,apy_old, apz_old), vfrac_old,
    //  //       AMREX_D_DECL(apx_new,apy_new, apz_new), vfrac_new,
    //  //       AMREX_D_DECL(fcx,fcy, fcz), ccent,
    //  //       d_bcrec_ptr, geom, dt, redistribution_type,
    //  //       defaults::srd_max_order, defaults::target_vol_fraction,
    //  //       amrex::Array4<amrex::Real const> {},
    //  //       vel_eb, bnorm, barea);
    // }

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
                              amrex::Real target_volfrac = 1.0);

    void StateRedistribute ( amrex::Box const& bx, int ncomp,
                             amrex::Array4<amrex::Real> const& dUdt_out,
                             amrex::Array4<amrex::Real> const& dUdt_in,
                             amrex::Array4<amrex::EBCellFlag const> const& flag,
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
                                amrex::Array4<amrex::EBCellFlag const> const& flag,
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
