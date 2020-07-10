#include <incflo_slopes_K.H>
#include <MOL.H>
#include <utility>

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir_or_ho (BCRec const* bcrec, int ncomp, int dir)
    {
        std::pair<bool,bool> r{false,false};
        for (int n = 0; n < ncomp; ++n) {
            r.first = r.first 
                 or (bcrec[n].lo(dir) == BCType::ext_dir)
                 or (bcrec[n].lo(dir) == BCType::hoextrap);
            r.second = r.second 
                 or (bcrec[n].hi(dir) == BCType::ext_dir)
                 or (bcrec[n].hi(dir) == BCType::hoextrap);
        }
        return r;
    }
}

#ifdef AMREX_USE_EB
void 
mol::compute_convective_fluxes_eb (int lev, Box const& bx, int ncomp,
                                   AMREX_D_DECL(Array4<Real> const& fx,
                                                Array4<Real> const& fy,
                                                Array4<Real> const& fz),
                                   Array4<Real const> const& q,
                                   AMREX_D_DECL(Array4<Real const> const& umac,
                                                Array4<Real const> const& vmac,
                                                Array4<Real const> const& wmac),
                                   BCRec const* h_bcrec,
                                   BCRec const* d_bcrec,
                                   Array4<EBCellFlag const> const& flag,
                                   AMREX_D_DECL(Array4<Real const> const& fcx,
                                                Array4<Real const> const& fcy,
                                                Array4<Real const> const& fcz),
                                   Array4<Real const> const& ccc,
                                   Vector<Geometry> geom)
{
    constexpr Real small_vel = 1.e-10;

    const Box& domain_box = geom[lev].Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);
#endif

    AMREX_D_TERM(Box const& xbx = amrex::surroundingNodes(bx,0);,
                 Box const& ybx = amrex::surroundingNodes(bx,1);,
                 Box const& zbx = amrex::surroundingNodes(bx,2););

    // ****************************************************************************
    // Predict to x-faces
    // ****************************************************************************

    // At an ext_dir or hoextrap boundary,
    //    the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir_or_ho(h_bcrec, ncomp, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo = extdir_lohi.first;
    bool has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_ilo >= xbx.smallEnd(0)-1) or
        (has_extdir_or_ho_hi and domain_ihi <= xbx.bigEnd(0)))
    {
        amrex::ParallelFor(xbx, ncomp,
        [d_bcrec,AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
        AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
        q,ccc,fcx,flag,umac,small_vel,fx]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           AMREX_D_TERM(bool extdir_or_ho_ilo = (d_bcrec[n].lo(0) == BCType::ext_dir) or
                                                (d_bcrec[n].lo(0) == BCType::hoextrap);,
                        bool extdir_or_ho_jlo = (d_bcrec[n].lo(1) == BCType::ext_dir) or
                                                (d_bcrec[n].lo(1) == BCType::hoextrap);,
                        bool extdir_or_ho_klo = (d_bcrec[n].lo(2) == BCType::ext_dir) or
                                                (d_bcrec[n].lo(2) == BCType::hoextrap););

           AMREX_D_TERM(bool extdir_or_ho_ihi = (d_bcrec[n].hi(0) == BCType::ext_dir) or
                                                (d_bcrec[n].hi(0) == BCType::hoextrap);,
                        bool extdir_or_ho_jhi = (d_bcrec[n].hi(1) == BCType::ext_dir) or
                                                (d_bcrec[n].hi(1) == BCType::hoextrap);,
                        bool extdir_or_ho_khi = (d_bcrec[n].hi(2) == BCType::ext_dir) or
                                                (d_bcrec[n].hi(2) == BCType::hoextrap););
           Real qs;

           if (flag(i,j,k).isConnected(-1,0,0)) 
           {
               if (i <= domain_ilo && (d_bcrec[n].lo(0) == BCType::ext_dir)) {
                   qs = q(domain_ilo-1,j,k,n);
               } else if (i >= domain_ihi+1 && (d_bcrec[n].hi(0) == BCType::ext_dir)) {
                   qs = q(domain_ihi+1,j,k,n);
               } else {

                   Real yf = fcx(i,j,k,0); // local (y,z) of centroid of z-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
                   Real zf = fcx(i,j,k,1);
#endif 

                   AMREX_D_TERM(Real xc = ccc(i,j,k,0);, // centroid of cell (i,j,k)
                                Real yc = ccc(i,j,k,1);,
                                Real zc = ccc(i,j,k,2););
 
                   AMREX_D_TERM(Real delta_x = 0.5 + xc;,
                                Real delta_y = yf  - yc;,
                                Real delta_z = zf  - zc;);
    
                   Real cc_qmax = amrex::max(q(i,j,k,n),q(i-1,j,k,n));
                   Real cc_qmin = amrex::min(q(i,j,k,n),q(i-1,j,k,n));
     
                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,n,q,ccc,flag,
                                              AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo), 
                                              AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi), 
                                              AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo), 
                                              AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

#if (AMREX_SPACEDIM == 3) 
                   Real qpls = q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                            + delta_y * slopes_eb_hi[1]
                                            + delta_z * slopes_eb_hi[2];
#else
                   Real qpls = q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                            + delta_y * slopes_eb_hi[1];
#endif
                   qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);
    
                   AMREX_D_TERM(xc = ccc(i-1,j,k,0);, // centroid of cell (i-1,j,k)
                                yc = ccc(i-1,j,k,1);,
                                zc = ccc(i-1,j,k,2););
    
                   AMREX_D_TERM(delta_x = 0.5 - xc;,
                                delta_y = yf  - yc;,
                                delta_z = zf  - zc;);

                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i-1,j,k,n,q,ccc,flag,
                                              AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                              AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                              AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                              AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

#if (AMREX_SPACEDIM == 3)    
                   Real qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                            + delta_y * slopes_eb_lo[1]
                                            + delta_z * slopes_eb_lo[2];
#else
                   Real qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                            + delta_y * slopes_eb_lo[1];
#endif
                   qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);

                    if (umac(i,j,k) > small_vel) {
                        qs = qmns;
                    } else if (umac(i,j,k) < -small_vel) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
               }

               fx(i,j,k,n) = umac(i,j,k) * qs;
   
           } else {
               fx(i,j,k,n) = 0.0;
           }
        });
    }
    else
    {
        amrex::ParallelFor(xbx, ncomp,
        [q,ccc,fcx,flag,umac,small_vel,fx]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           Real qs;
           if (flag(i,j,k).isConnected(-1,0,0)) 
           {
               Real yf = fcx(i,j,k,0); // local (y,z) of centroid of z-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcx(i,j,k,1);
#endif

               AMREX_D_TERM(Real xc = ccc(i,j,k,0);, // centroid of cell (i,j,k)
                            Real yc = ccc(i,j,k,1);,
                            Real zc = ccc(i,j,k,2););

               AMREX_D_TERM(Real delta_x = 0.5 + xc;,
                            Real delta_y = yf  - yc;,
                            Real delta_z = zf  - zc;);

               Real cc_qmax = amrex::max(q(i,j,k,n),q(i-1,j,k,n));
               Real cc_qmin = amrex::min(q(i,j,k,n),q(i-1,j,k,n));

               // Compute slopes of component "n" of q
               const auto& slopes_eb_hi = incflo_slopes_eb(i,j,k,n,q,ccc,flag);

#if (AMREX_SPACEDIM == 3)
               Real qpls = q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                        + delta_y * slopes_eb_hi[1]
                                        + delta_z * slopes_eb_hi[2];
#else
               Real qpls = q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                        + delta_y * slopes_eb_hi[1];
#endif
               qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);

               AMREX_D_TERM(xc = ccc(i-1,j,k,0);, // centroid of cell (i-1,j,k)
                            yc = ccc(i-1,j,k,1);,
                            zc = ccc(i-1,j,k,2););

               AMREX_D_TERM(delta_x = 0.5 - xc;,
                            delta_y = yf  - yc;,
                            delta_z = zf  - zc;);

               // Compute slopes of component "n" of q
               const auto& slopes_eb_lo = incflo_slopes_eb(i-1,j,k,n,q,ccc,flag);

#if (AMREX_SPACEDIM == 3)
               Real qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                        + delta_y * slopes_eb_lo[1]
                                        + delta_z * slopes_eb_lo[2];
#else
               Real qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                        + delta_y * slopes_eb_lo[1];
#endif
               qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);

               if (umac(i,j,k) > small_vel) {
                    qs = qmns;
                } else if (umac(i,j,k) < -small_vel) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }

                fx(i,j,k,n) = umac(i,j,k) * qs;

           } else {
               fx(i,j,k,n) = 0.0;
           }
        });
    }

    // ****************************************************************************
    // Predict to y-faces
    // ****************************************************************************

    // At an ext_dir or hoextrap boundary,
    //    the boundary value is on the face, not cell center.
    extdir_lohi = has_extdir_or_ho(h_bcrec, ncomp,  static_cast<int>(Direction::y));
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_jlo >= ybx.smallEnd(1)-1) or
        (has_extdir_or_ho_hi and domain_jhi <= ybx.bigEnd(1)))
    {
        amrex::ParallelFor(ybx, ncomp,
        [d_bcrec,AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
         AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
         q,ccc,fcy,flag,vmac,small_vel,fy]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qs;
            if (flag(i,j,k).isConnected(0,-1,0)) 
            {
                AMREX_D_TERM(bool extdir_or_ho_ilo = (d_bcrec[n].lo(0) == BCType::ext_dir) or
                                                     (d_bcrec[n].lo(0) == BCType::hoextrap);,
                             bool extdir_or_ho_jlo = (d_bcrec[n].lo(1) == BCType::ext_dir) or
                                                     (d_bcrec[n].lo(1) == BCType::hoextrap);,
                             bool extdir_or_ho_klo = (d_bcrec[n].lo(2) == BCType::ext_dir) or
                                                     (d_bcrec[n].lo(2) == BCType::hoextrap););
                AMREX_D_TERM(bool extdir_or_ho_ihi = (d_bcrec[n].hi(0) == BCType::ext_dir) or
                                                     (d_bcrec[n].hi(0) == BCType::hoextrap);,
                             bool extdir_or_ho_jhi = (d_bcrec[n].hi(1) == BCType::ext_dir) or
                                                     (d_bcrec[n].hi(1) == BCType::hoextrap);,
                             bool extdir_or_ho_khi = (d_bcrec[n].hi(2) == BCType::ext_dir) or
                                                     (d_bcrec[n].hi(2) == BCType::hoextrap););

                if (j <= domain_jlo && (d_bcrec[n].lo(1) == BCType::ext_dir)) {
                    qs = q(i,domain_jlo-1,k,n);
                } else if (j >= domain_jhi+1 && (d_bcrec[n].hi(1) == BCType::ext_dir)) {
                    qs = q(i,domain_jhi+1,k,n);
                } else {

                   Real xf = fcy(i,j,k,0); // local (x,z) of centroid of z-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
                   Real zf = fcy(i,j,k,1);
#endif

                   AMREX_D_TERM(Real xc = ccc(i,j,k,0);, // centroid of cell (i,j,k)
                                Real yc = ccc(i,j,k,1);,
                                Real zc = ccc(i,j,k,2););
 
                   AMREX_D_TERM(Real delta_x = xf  - xc;,
                                Real delta_y = 0.5 + yc;,
                                Real delta_z = zf  - zc;);
    
                   Real cc_qmax = amrex::max(q(i,j,k,n),q(i,j-1,k,n));
                   Real cc_qmin = amrex::min(q(i,j,k,n),q(i,j-1,k,n));
     
                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,n,q,ccc,flag,
                                              AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                              AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                              AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                              AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));
 
#if (AMREX_SPACEDIM == 3)
                   Real qpls = q(i,j  ,k,n) + delta_x * slopes_eb_hi[0]
                                            - delta_y * slopes_eb_hi[1]
                                            + delta_z * slopes_eb_hi[2];
#else
                   Real qpls = q(i,j  ,k,n) + delta_x * slopes_eb_hi[0]
                                            - delta_y * slopes_eb_hi[1];
#endif
                   qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);
    
                   AMREX_D_TERM(xc = ccc(i,j-1,k,0);, // centroid of cell (i-1,j,k)
                                yc = ccc(i,j-1,k,1);,
                                zc = ccc(i,j-1,k,2););
    
                   AMREX_D_TERM(delta_x = xf  - xc;,
                                delta_y = 0.5 - yc;,
                                delta_z = zf  - zc;);

                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i,j-1,k,n,q,ccc,flag,
                                              AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                              AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                              AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                              AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

#if (AMREX_SPACEDIM == 3)    
                   Real qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                            + delta_y * slopes_eb_lo[1]
                                            + delta_z * slopes_eb_lo[2];
#else
                   Real qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                            + delta_y * slopes_eb_lo[1];
#endif
                   qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);

                    if (vmac(i,j,k) > small_vel) {
                        qs = qmns;
                    } else if (vmac(i,j,k) < -small_vel) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
                }

                fy(i,j,k,n) = vmac(i,j,k) * qs;

           } else {
                fy(i,j,k,n) = 0.0;
           }
        });
    }
    else
    {
        amrex::ParallelFor(ybx, ncomp,
        [q,ccc,fcy,flag,vmac,small_vel,fy]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qs;
            if (flag(i,j,k).isConnected(0,-1,0)) {

               Real xf = fcy(i,j,k,0); // local (x,z) of centroid of z-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcy(i,j,k,1);
#endif

               AMREX_D_TERM(Real xc = ccc(i,j,k,0);, // centroid of cell (i,j,k)
                            Real yc = ccc(i,j,k,1);,
                            Real zc = ccc(i,j,k,2););

               AMREX_D_TERM(Real delta_x = xf  - xc;,
                            Real delta_y = 0.5 + yc;,
                            Real delta_z = zf  - zc;);

               Real cc_qmax = amrex::max(q(i,j,k,n),q(i,j-1,k,n));
               Real cc_qmin = amrex::min(q(i,j,k,n),q(i,j-1,k,n));

               // Compute slopes of component "n" of q
               const auto& slopes_eb_hi = incflo_slopes_eb(i,j,k,n,q,ccc,flag);

#if (AMREX_SPACEDIM == 3)
               Real qpls = q(i,j  ,k,n) + delta_x * slopes_eb_hi[0]
                                        - delta_y * slopes_eb_hi[1]
                                        + delta_z * slopes_eb_hi[2];
#else
               Real qpls = q(i,j  ,k,n) + delta_x * slopes_eb_hi[0]
                                        - delta_y * slopes_eb_hi[1];
#endif
               qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);

               AMREX_D_TERM(xc = ccc(i,j-1,k,0);, // centroid of cell (i-1,j,k)
                            yc = ccc(i,j-1,k,1);,
                            zc = ccc(i,j-1,k,2););

               AMREX_D_TERM(delta_x = xf  - xc;,
                            delta_y = 0.5 - yc;,
                            delta_z = zf  - zc;);

               // Compute slopes of component "n" of q
               const auto& slopes_eb_lo = incflo_slopes_eb(i,j-1,k,n,q,ccc,flag);

#if (AMREX_SPACEDIM == 3)
               Real qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                        + delta_y * slopes_eb_lo[1]
                                        + delta_z * slopes_eb_lo[2];
#else
               Real qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                        + delta_y * slopes_eb_lo[1];
#endif
               qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);

               if (vmac(i,j,k) > small_vel) {
                   qs = qmns;
               } else if (vmac(i,j,k) < -small_vel) {
                   qs = qpls;
               } else {
                   qs = 0.5*(qmns+qpls);
               }

               fy(i,j,k,n) = vmac(i,j,k) * qs;

           } else {
               fy(i,j,k,n) = 0.0;
           }
        });
    }

    // ****************************************************************************
    // Predict to z-faces
    // ****************************************************************************
#if (AMREX_SPACEDIM == 3)
    // At an ext_dir or hoextrap boundary,
    //    the boundary value is on the face, not cell center.
    extdir_lohi = has_extdir_or_ho(h_bcrec, ncomp, static_cast<int>(Direction::z));
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_klo >= zbx.smallEnd(2)-1) or
        (has_extdir_or_ho_hi and domain_khi <= zbx.bigEnd(2)))
    {
        amrex::ParallelFor(zbx, ncomp,
        [d_bcrec,domain_ilo,domain_jlo,domain_klo,domain_ihi,domain_jhi,domain_khi,
         q,ccc,fcz,flag,wmac,small_vel,fz]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,0,-1)) {

                bool extdir_or_ho_ilo = (d_bcrec[n].lo(0) == BCType::ext_dir) or
                                        (d_bcrec[n].lo(0) == BCType::hoextrap);
                bool extdir_or_ho_ihi = (d_bcrec[n].hi(0) == BCType::ext_dir) or
                                        (d_bcrec[n].hi(0) == BCType::hoextrap);
                bool extdir_or_ho_jlo = (d_bcrec[n].lo(1) == BCType::ext_dir) or
                                        (d_bcrec[n].lo(1) == BCType::hoextrap);
                bool extdir_or_ho_jhi = (d_bcrec[n].hi(1) == BCType::ext_dir) or
                                        (d_bcrec[n].hi(1) == BCType::hoextrap);
                bool extdir_or_ho_klo = (d_bcrec[n].lo(2) == BCType::ext_dir) or
                                        (d_bcrec[n].lo(2) == BCType::hoextrap);
                bool extdir_or_ho_khi = (d_bcrec[n].hi(2) == BCType::ext_dir) or
                                        (d_bcrec[n].hi(2) == BCType::hoextrap);

                Real qs;
                if (k <= domain_klo && (d_bcrec[n].lo(2) == BCType::ext_dir)) {
                    qs = q(i,j,domain_klo-1,n);
                } else if (k >= domain_khi+1 && (d_bcrec[n].hi(2) == BCType::ext_dir)) {
                    qs = q(i,j,domain_khi+1,n);
                } else {

                    Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
                    Real yf = fcz(i,j,k,1);
 
                    Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
                    Real yc = ccc(i,j,k,1);
                    Real zc = ccc(i,j,k,2);
 
                    Real delta_x = xf  - xc;
                    Real delta_y = yf  - yc;
                    Real delta_z = 0.5 + zc;
     
                    Real cc_qmax = amrex::max(q(i,j,k,n),q(i,j,k-1,n));
                    Real cc_qmin = amrex::min(q(i,j,k,n),q(i,j,k-1,n));
     
                    // Compute slopes of component "n" of q
                    const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,n,q,ccc,flag,
                                               AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                               AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                               AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                               AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));
 
                    Real qpls = q(i,j,k  ,n) + delta_x * slopes_eb_hi[0]
                                             + delta_y * slopes_eb_hi[1]
                                             - delta_z * slopes_eb_hi[2];
     
                    qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);
     
                    xc = ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
                    yc = ccc(i,j,k-1,1);
                    zc = ccc(i,j,k-1,2);
     
                    delta_x = xf  - xc;
                    delta_y = yf  - yc;
                    delta_z = 0.5 - zc;

                    // Compute slopes of component "n" of q
                    const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i,j,k-1,n,q,ccc,flag,
                                               AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                               AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                               AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                               AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));
    
                    Real qmns = q(i,j,k-1,n) + delta_x * slopes_eb_lo[0]
                                             + delta_y * slopes_eb_lo[1]
                                             + delta_z * slopes_eb_lo[2];
    
                    qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);

                    if (wmac(i,j,k) > small_vel) {
                        qs = qmns;
                    } else if (wmac(i,j,k) < -small_vel) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
                }

                fz(i,j,k,n) = wmac(i,j,k) * qs;

           } else {
                fz(i,j,k,n) = 0.0;
           }
        });
    }
    else
    {
        amrex::ParallelFor(zbx, ncomp,
        [q,ccc,fcz,flag,wmac,small_vel,fz]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,0,-1)) {
                Real qs;

                Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
                Real yf = fcz(i,j,k,1);
 
                Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
                Real yc = ccc(i,j,k,1);
                Real zc = ccc(i,j,k,2);
 
                Real delta_x = xf  - xc;
                Real delta_y = yf  - yc;
                Real delta_z = 0.5 + zc;
     
                Real cc_qmax = amrex::max(q(i,j,k,n),q(i,j,k-1,n));
                Real cc_qmin = amrex::min(q(i,j,k,n),q(i,j,k-1,n));
     
                // Compute slopes of component "n" of q
                const auto& slopes_eb_hi = incflo_slopes_eb(i,j,k,n,q,ccc,flag);
 
                Real qpls = q(i,j,k  ,n) + delta_x * slopes_eb_hi[0]
                                         + delta_y * slopes_eb_hi[1]
                                         - delta_z * slopes_eb_hi[2];
 
                qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);
 
                xc = ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
                yc = ccc(i,j,k-1,1);
                zc = ccc(i,j,k-1,2);
 
                delta_x = xf  - xc;
                delta_y = yf  - yc;
                delta_z = 0.5 - zc;

                // Compute slopes of component "n" of q
                const auto& slopes_eb_lo = incflo_slopes_eb(i,j,k-1,n,q,ccc,flag);
    
                Real qmns = q(i,j,k-1,n) + delta_x * slopes_eb_lo[0]
                                         + delta_y * slopes_eb_lo[1]
                                         + delta_z * slopes_eb_lo[2];

                qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);

                if (wmac(i,j,k) > small_vel) {
                    qs = qmns;
                } else if (wmac(i,j,k) < -small_vel) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }

                fz(i,j,k,n) = wmac(i,j,k) * qs;

           } else {
                fz(i,j,k,n) = 0.0;
           }
        });
    }
#endif
}
#endif
