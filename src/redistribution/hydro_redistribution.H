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

/**new
 * Placeholder description of Redistribution namespace.
 *
 */

constexpr amrex::Real target_vf = 1.0;

namespace Redistribution {

    void Apply ( amrex::Box const& bx, int ncomp,
                 amrex::Array4<amrex::Real>       const& dUdt_out,
                 amrex::Array4<amrex::Real>       const& dUdt_in,
                 amrex::Array4<amrex::Real const> const& U_in,
                 amrex::Array4<amrex::Real> const& scratch,
                 amrex::Array4<amrex::EBCellFlag const> const& flag_old,
		 amrex::Array4<amrex::EBCellFlag const> const& flag_new,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx_new,
                              amrex::Array4<amrex::Real const> const& apy_new,
                              amrex::Array4<amrex::Real const> const& apz_new),
		 amrex::Array4<amrex::Real const> const& vfrac_old,
                 amrex::Array4<amrex::Real const> const& vfrac_new,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                              amrex::Array4<amrex::Real const> const& fcy,
                              amrex::Array4<amrex::Real const> const& fcz),
                 amrex::Array4<amrex::Real const> const& ccent,
                 amrex::BCRec  const* d_bcrec_ptr,
                 amrex::Geometry const& geom,
                 amrex::Real dt, std::string redistribution_type,
                 const int srd_max_order = 2,
                 amrex::Real target_volfrac = target_vf,
                 amrex::Array4<amrex::Real const> const& update_scale={});

    void ApplyToInitialData ( amrex::Box const& bx, int ncomp,
                              amrex::Array4<amrex::Real                  > const& U_out,
                              amrex::Array4<amrex::Real                  > const& U_in,
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
                              amrex::Geometry& geom, std::string redistribution_type,
                              const int srd_max_order = 2,
                              amrex::Real target_volfrac = target_vf);

    void StateRedistribute ( amrex::Box const& bx, int ncomp,
                             amrex::Array4<amrex::Real> const& dUdt_out,
                             amrex::Array4<amrex::Real> const& dUdt_in,
                             amrex::Array4<amrex::EBCellFlag const> const& flag_old,
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
                             const int max_order = 2);

    void MakeITracker ( amrex::Box const& bx,
			amrex::Array4<amrex::EBCellFlag const> const& flag_old,
			amrex::Array4<amrex::EBCellFlag const> const& flag_new,
			AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx_new,
                                     amrex::Array4<amrex::Real const> const& apy_new,
                                     amrex::Array4<amrex::Real const> const& apz_new),
                        amrex::Array4<amrex::Real const> const& vfrac_old,
			amrex::Array4<amrex::Real const> const& vfrac_new,
                        amrex::Array4<int> const& itracker,
                        amrex::Geometry const& geom,
                        amrex::Real target_volfrac);

    void MakeStateRedistUtils ( amrex::Box const& bx,
                                amrex::Array4<amrex::EBCellFlag const> const& flag_old,
				amrex::Array4<amrex::EBCellFlag const> const& flag_new,
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
