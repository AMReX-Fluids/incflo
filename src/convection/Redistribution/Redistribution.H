//
//  Redistribution.H
//

#ifndef Redistribution_H
#define Redistribution_H

#ifdef AMREX_USE_EB

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiCutFab.H>

namespace redistribution {

    void redistribute_eb (amrex::Box const& bx, int ncomp,
                          amrex::Array4<amrex::Real>       const& dUdt_out,
                          amrex::Array4<amrex::Real>       const& dUdt_in,
                          amrex::Array4<amrex::Real const> const& U_in,
                          amrex::Array4<amrex::Real> const& scratch,
                          AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                       amrex::Array4<amrex::Real const> const& vmac,
                                       amrex::Array4<amrex::Real const> const& wmac),
                          amrex::Array4<amrex::EBCellFlag const> const& flag,
                          AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                       amrex::Array4<amrex::Real const> const& apy,
                                       amrex::Array4<amrex::Real const> const& apz),
                          amrex::Array4<amrex::Real const> const& vfrac,
                          AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                       amrex::Array4<amrex::Real const> const& fcy,
                                       amrex::Array4<amrex::Real const> const& fcz),
                          amrex::Array4<amrex::Real const> const& ccent,
                          amrex::Geometry& geom, amrex::Real dt, std::string redistribution_type);

    void redistribute_initial_data (amrex::Box const& bx, int ncomp,
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
                          amrex::Geometry& geom, std::string redistribution_type);

    void flux_redistribute_eb ( amrex::Box const& bx, int ncomp,
                                amrex::Array4<amrex::Real>       const& dUdt_out,
                                amrex::Array4<amrex::Real const> const& dUdt_in,
                                amrex::Array4<amrex::Real> const& scratch,
                                amrex::Array4<amrex::EBCellFlag const> const& flag,
                                amrex::Array4<amrex::Real const> const& vfrac,
                                amrex::Geometry& geom);

    void state_redistribute( amrex::Box const& bx, int ncomp,
                             amrex::Array4<amrex::Real> const& dUdt_out,
                             amrex::Array4<amrex::Real> const& dUdt_in,
                             amrex::Array4<amrex::EBCellFlag const> const& flag,
                             AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                          amrex::Array4<amrex::Real const> const& apy,
                                          amrex::Array4<amrex::Real const> const& apz),
                             amrex::Array4<amrex::Real const> const& vfrac,
                             AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                          amrex::Array4<amrex::Real const> const& fcy,
                                          amrex::Array4<amrex::Real const> const& fcz),
                             amrex::Array4<amrex::Real const> const& ccent,
                             amrex::Geometry& geom);

    void make_itracker             ( amrex::Box const& bx,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                                  amrex::Array4<amrex::Real const> const& apy,
                                                  amrex::Array4<amrex::Real const> const& apz),
                                     amrex::Array4<amrex::Real const> const& vfrac,
                                     amrex::Array4<int> const& itracker,
                                     amrex::Geometry& geom);

    void merge_redistribute( amrex::Box const& bx, int ncomp,
                             amrex::Array4<amrex::Real> const& dUdt_out,
                             amrex::Array4<amrex::Real> const& dUdt_in,
                             AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                          amrex::Array4<amrex::Real const> const& apy,
                                         amrex::Array4<amrex::Real const> const& apz),
                            amrex::Array4<amrex::Real const> const& vfrac,
                            amrex::Array4<int> const& itracker,
                            amrex::Geometry& geom);

} // namespace redistribution

#endif

#endif /* Redistribution_H */
