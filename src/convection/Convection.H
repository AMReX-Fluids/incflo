//
//  Convection.H
//

#ifndef Convection_H
#define Convection_H

#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>

#ifdef AMREX_USE_EB
#include <AMReX_MultiCutFab.H>
#endif

namespace convection {

#ifdef AMREX_USE_EB

    void redistribute_convective_term (amrex::Box const& bx, amrex::MFIter const& mfi,
                                       amrex::Array4<amrex::Real const> const& vel,       // velocity
                                       amrex::Array4<amrex::Real const> const& rho,       // density
                                       amrex::Array4<amrex::Real const> const& rhotrac,   // tracer
                                       amrex::Array4<amrex::Real> const& dvdt_tmp,  // initial velocity update
                                       amrex::Array4<amrex::Real> const& drdt_tmp,  // initial density update
                                       amrex::Array4<amrex::Real> const& dtdt_tmp,  // initial tracer update
                                       amrex::Array4<amrex::Real> const& dvdt,      // final velocity update
                                       amrex::Array4<amrex::Real> const& drdt,      // final density update
                                       amrex::Array4<amrex::Real> const& dtdt,      // final tracer update
                                       std::string l_redistribution_type, 
                                       bool l_constant_density, bool l_advect_tracer, int l_ntrac,
                                       amrex::EBFArrayBoxFactory const* ebfact,
                                       amrex::Geometry& geom, amrex::Real l_dt);
#endif
} // namespace convection

#endif /* Convection_H */
