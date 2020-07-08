//
//  MOL.H
//

#ifndef MOL_H
#define MOL_H

#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>

#ifdef AMREX_USE_EB
#include <AMReX_MultiCutFab.H>
#endif


namespace mol {

    void predict_vels_on_faces (int lev,  
                                AMREX_D_DECL(amrex::MultiFab& u_mac,
                                             amrex::MultiFab& v_mac, 
                                             amrex::MultiFab& w_mac),
                                amrex::MultiFab const& vel,
                                amrex::Vector<amrex::BCRec> const& h_bcrec,
                                              amrex::BCRec  const* d_bcrec,
#ifdef AMREX_USE_EB
                                amrex::EBFArrayBoxFactory const* ebfact,
#endif
                                amrex::Vector<amrex::Geometry> geom);

    void predict_vels_on_faces (int lev,
                                AMREX_D_DECL(amrex::Box const& ubx,
                                             amrex::Box const& vbx,
                                             amrex::Box const& wbx),
                                AMREX_D_DECL(amrex::Array4<amrex::Real> const& u,
                                             amrex::Array4<amrex::Real> const& v,
                                             amrex::Array4<amrex::Real> const& w),
                                amrex::Array4<amrex::Real const> const& vcc,
                                amrex::Vector<amrex::BCRec> const& h_bcrec,
                                              amrex::BCRec  const* d_bcrec,
                                amrex::Vector<amrex::Geometry> geom);

    void compute_convective_fluxes (int lev, amrex::Box const& bx, int ncomp,
                                    AMREX_D_DECL(amrex::Array4<amrex::Real> const& fx,
                                                 amrex::Array4<amrex::Real> const& fy,
                                                 amrex::Array4<amrex::Real> const& fz),
                                    amrex::Array4<amrex::Real const> const& q,
                                    AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                                 amrex::Array4<amrex::Real const> const& vmac,
                                                 amrex::Array4<amrex::Real const> const& wmac),
                                    amrex::BCRec const* h_bcrec,
                                    amrex::BCRec const* d_bcrec,
                                    amrex::Vector<amrex::Geometry> geom);

    void compute_convective_rate (int lev, amrex::Box const& bx, int ncomp,
                                  amrex::Array4<amrex::Real> const& dUdt,
                                  AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fx,
                                               amrex::Array4<amrex::Real const> const& fy,
                                               amrex::Array4<amrex::Real const> const& fz),
                                  amrex::Vector<amrex::Geometry> geom);

#ifdef AMREX_USE_EB
    void predict_vels_on_faces_eb (int lev,
                                   amrex::Box const& ccbx,
                                   AMREX_D_DECL(amrex::Box const& ubx,
                                                amrex::Box const& vbx,
                                                amrex::Box const& wbx),
                                   AMREX_D_DECL(amrex::Array4<amrex::Real> const& u,
                                                amrex::Array4<amrex::Real> const& v,
                                                amrex::Array4<amrex::Real> const& w),
                                   amrex::Array4<amrex::Real const> const& vcc,
                                   amrex::Array4<amrex::EBCellFlag const> const& flag,
                                   AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                                amrex::Array4<amrex::Real const> const& fcy,
                                                amrex::Array4<amrex::Real const> const& fcz), 
                                   amrex::Array4<amrex::Real const> const& ccc,
                                   amrex::Vector<amrex::BCRec> const& h_bcrec,
                                                 amrex::BCRec  const* d_bcrec,
                                   amrex::Vector<amrex::Geometry> geom);

    void compute_convective_fluxes_eb (int lev, amrex::Box const& bx, int ncomp,
                                       AMREX_D_DECL(amrex::Array4<amrex::Real> const& fx,
                                                    amrex::Array4<amrex::Real> const& fy,
                                                    amrex::Array4<amrex::Real> const& fz),
                                       amrex::Array4<amrex::Real const> const& q,
                                       AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                                    amrex::Array4<amrex::Real const> const& vmac,
                                                    amrex::Array4<amrex::Real const> const& wmac),
                                       amrex::BCRec const* h_bcrec,
                                       amrex::BCRec const* d_bcrec,
                                       amrex::Array4<amrex::EBCellFlag const> const& flag,
                                       AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                                    amrex::Array4<amrex::Real const> const& fcy,
                                                    amrex::Array4<amrex::Real const> const& fcz),
                                       amrex::Array4<amrex::Real const> const& ccc,
                                       amrex::Vector<amrex::Geometry> geom);

    void compute_convective_rate_eb (int lev, amrex::Box const& bx, int ncomp,
                                     amrex::Array4<amrex::Real> const& dUdt,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fx,
                                                  amrex::Array4<amrex::Real const> const& fy,
                                                  amrex::Array4<amrex::Real const> const& fz),
                                     amrex::Array4<amrex::EBCellFlag const> const& flag,
                                     amrex::Array4<amrex::Real const> const& vfrac,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                                  amrex::Array4<amrex::Real const> const& apy,
                                                  amrex::Array4<amrex::Real const> const& apz),
                                     amrex::Vector<amrex::Geometry> geom);
#endif
} // namespace mol

#endif /* MOL_H */
