#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <incflo.H>

using namespace amrex;

// Make the iMultiFab needed by the Nodal Projection to support mixed BCs
iMultiFab
incflo::make_nodalBC_mask(int lev)
{
    // We use the solver's overset mask to implement mixed BCs. The mask is nodal, and
    // so exists on the boundary face, where 0 indicates a Dirichlet boundary condition
    // (i,e, no solve needed for that point because we know the solution there already).
    // Because it's an overset mask and not a BC, we must fill the entire MF.
    iMultiFab new_mask(amrex::convert(grids[lev], IntVect::TheNodeVector()), dmap[lev], 1, 0);
    new_mask = 1;
    int inflow = 1;
    int outflow = 0;

    Geometry const& gm = Geom(lev);
    Box const& domain = gm.Domain();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        Orientation olo(dir,Orientation::low);
        Orientation ohi(dir,Orientation::high);
        if (m_bc_type[olo] == BC::mixed || m_bc_type[ohi] == BC::mixed) {
            Box dlo = (m_bc_type[olo] == BC::mixed) ? surroundingNodes(bdryLo(domain,dir)) : Box();
            Box dhi = (m_bc_type[ohi] == BC::mixed) ? surroundingNodes(bdryHi(domain,dir)) : Box();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(new_mask); mfi.isValid(); ++mfi) {
                Box blo = mfi.validbox() & dlo;
                Box bhi = mfi.validbox() & dhi;
                Array4<int> const& mask_arr = new_mask.array(mfi);
                if (blo.ok()) {
                    prob_set_BC_MF(olo, blo, mask_arr, lev, inflow, outflow, "projection");
                }
                if (bhi.ok()) {
                    prob_set_BC_MF(ohi, bhi, mask_arr, lev, inflow, outflow, "projection");
                }
            }
        }
    }

    return new_mask;
}

// Make a position dependent BC MultiFab for the advection routines to support
// mixed BCs.
std::unique_ptr<iMultiFab>
incflo::make_BC_MF(int lev, amrex::Gpu::DeviceVector<amrex::BCRec> const& bcs,
                   std::string const& field)
{
    auto ncomp = static_cast<int>(bcs.size());
    // The advection routines expect that the BC type is stored in the first ghost
    // cell of a cell-centered MF.
    std::unique_ptr<iMultiFab> BC_MF(new iMultiFab(grids[lev], dmap[lev], ncomp, 1));
    // initialize to 0 == BCType::int_dir, not sure this is needed really
    *BC_MF = 0;
    // For advection, we use FOEXTRAP for outflow regions per incflo::init_bcs()
    int inflow = BCType::ext_dir;
    int outflow = BCType::foextrap;

    Box const& domain = geom[lev].Domain();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        Orientation olo(dir,Orientation::low);
        Orientation ohi(dir,Orientation::high);

        Box dlo = adjCellLo(domain,dir);
        Box dhi = adjCellHi(domain,dir);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*BC_MF); mfi.isValid(); ++mfi) {
            Box blo = mfi.growntilebox() & dlo;
            Box bhi = mfi.growntilebox() & dhi;
            Array4<int> const& bc_arr = BC_MF->array(mfi);
            const auto bc_ptr = bcs.data();
            if (m_bc_type[olo] == BC::mixed) {
                prob_set_BC_MF(olo, blo, bc_arr, lev, inflow, outflow, field);
            } else {
                amrex::ParallelFor(blo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int n = 0; n < ncomp; n++) {
                        bc_arr(i,j,k,n) = bc_ptr[n].lo(dir);
                    }
                });
            }
            if (m_bc_type[ohi] == BC::mixed) {
                prob_set_BC_MF(ohi, bhi, bc_arr, lev, inflow, outflow, field);
            } else {
                amrex::ParallelFor(bhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int n = 0; n < ncomp; n++) {
                        bc_arr(i,j,k,n) = bc_ptr[n].hi(dir);
                    }
                });
            }
        }
    }

    return BC_MF;
}

// Create the MutliFabs needed by the linear solver to define Robin BCs, which incflo
// uses to support mixed BCs in the MAC Projection and diffusion.
Vector<MultiFab>
incflo::make_robinBC_MFs(int lev, MultiFab* state)
{
    // MLMG defines the Robin BC with
    //     a u + b du/dn = f
    // The values are stored in the first ghost cell, although the BC is considered on
    // the face. Only the ghost cells of robin BC arrays are used in MLMG.
    int nghost = 1;

    Vector<MultiFab> robin(3);
    for (int n = 0; n < 3; n++) {
        robin[n].define(grids[lev], dmap[lev], 1, nghost);
    }

    MultiFab& robin_a = robin[0];
    MultiFab& robin_b = robin[1];
    MultiFab& robin_f = robin[2];

    Box const& domain = Geom(lev).Domain();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        Orientation olo(dir,Orientation::low);
        Orientation ohi(dir,Orientation::high);

        // Only need to fill Robin BC sides, MLMG will check for Robin BC first
        if (m_bc_type[olo] == BC::mixed || m_bc_type[ohi] == BC::mixed) {
            Box dlo = (m_bc_type[olo] == BC::mixed) ? adjCellLo(domain,dir) : Box();
            Box dhi = (m_bc_type[ohi] == BC::mixed) ? adjCellHi(domain,dir) : Box();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            // FIXME - do we want tiling here...
            for (MFIter mfi(robin_a,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& gbx = amrex::grow(mfi.validbox(),nghost);
                Box blo = gbx & dlo;
                Box bhi = gbx & dhi;
                Array4<Real> const& a_arr     = robin_a.array(mfi);
                Array4<Real> const& b_arr     = robin_b.array(mfi);
                Array4<Real> const& f_arr     = robin_f.array(mfi);

                if (state) {
                    // state has filled boundry, so can use its ghost cells for the BC values
                    Array4<Real const> const& bcv = state->const_array(mfi);
                    // Diffusion, solving for velocity, density, tracer
                    //
                    // Robin BC:   a u + b du/dn = f  -- inflow,  Dirichlet a=1, b=0, f=bcv
                    //                                -- outflow, Neumann   a=0, b=1, f=0
                    if (blo.ok()) {
                        prob_set_diffusion_robinBCs(olo, blo, a_arr, b_arr, f_arr, bcv, lev);
                    }
                    if (bhi.ok()) {
                        prob_set_diffusion_robinBCs(ohi, bhi, a_arr, b_arr, f_arr, bcv, lev);
                    }
                } else {
                    // MAC, solving for phi (~pressure)
                    //
                    // Robin BC:   a u + b du/dn = f  -- inflow,  Neumann   a=0, b=1, f=0
                    //                                -- outflow, Dirichlet a=1, b=0, f=0
                    if (blo.ok()) {
                        prob_set_MAC_robinBCs(olo, blo, a_arr, b_arr, f_arr, lev);
                    }
                    if (bhi.ok()) {
                        prob_set_MAC_robinBCs(ohi, bhi, a_arr, b_arr, f_arr, lev);
                    }
                }
            }
        }
    }

    return robin;
}

#ifdef AMREX_USE_EB
void
incflo::set_eb_velocity (int lev, amrex::Real /*time*/, MultiFab& eb_vel, int nghost)
{
    Geometry const& gm = Geom(lev);
    eb_vel.setVal(0.);

    const auto& factory =
       dynamic_cast<EBFArrayBoxFactory const&>(eb_vel.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(eb_vel, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       const Box& bx = mfi.tilebox();
       const auto& flagfab      = factory.getMultiEBCellFlagFab()[mfi];

       if (flagfab.getType(bx) == FabType::singlevalued) {
          const auto& flags_arr    = flagfab.const_array();
          const auto& eb_vel_arr   = eb_vel[mfi].array();
          const auto& norm_arr     = factory.getBndryNormal()[mfi].const_array();

          bool has_normal = m_eb_flow.has_normal;
          GpuArray<Real, AMREX_SPACEDIM> normal{0.};
          if (has_normal) {
             AMREX_D_TERM(
             normal[0] = m_eb_flow.normal[0];,
             normal[1] = m_eb_flow.normal[1];,
             normal[2] = m_eb_flow.normal[2]);
          }
          Real pad = std::numeric_limits<float>::epsilon();
          Real normal_tol = m_eb_flow.normal_tol;
          Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
          Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

          bool has_comps = false;
          Real eb_vel_mag(0.0);
          GpuArray<amrex::Real,3> eb_vel_comps{0.};
          // Flow is specified as a velocity magnitude
          if ( m_eb_flow.is_mag ) {
             eb_vel_mag = m_eb_flow.vel_mag;
          } else {
             has_comps = true;
             const auto& vels = m_eb_flow.velocity;
             AMREX_D_TERM(eb_vel_comps[0] = vels[0];,
                          eb_vel_comps[1] = vels[1];,
                          eb_vel_comps[2] = vels[2]);
          }

          ParallelFor(bx, [flags_arr,eb_vel_arr,norm_arr,has_comps,has_normal,normal,
                 norm_tol_lo, norm_tol_hi,eb_vel_mag,eb_vel_comps]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
             if (flags_arr(i,j,k).isSingleValued()) {
                Real mask = Real(1.0);

                if(has_normal) {
#if (AMREX_SPACEDIM==3)
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1]
                                 + norm_arr(i,j,k,2)*normal[2];
#else
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1];
#endif

                  mask = ((norm_tol_lo <= dotprod) &&
                          (dotprod <= norm_tol_hi)) ? Real(1.0) : Real(0.0);
                }

                if (has_comps) {
                   AMREX_D_TERM(eb_vel_arr(i,j,k,0) = mask*eb_vel_comps[0];,
                                eb_vel_arr(i,j,k,1) = mask*eb_vel_comps[1];,
                                eb_vel_arr(i,j,k,2) = mask*eb_vel_comps[2]);
                } else {
                   // The EB normal points out of the domain so we need to flip the
                   // when using it to convert magnitude to velocity components so
                   // the resulting vector points into the domain.
                   AMREX_D_TERM(eb_vel_arr(i,j,k,0) = -mask*norm_arr(i,j,k,0)*eb_vel_mag;,
                                eb_vel_arr(i,j,k,1) = -mask*norm_arr(i,j,k,1)*eb_vel_mag;,
                                eb_vel_arr(i,j,k,2) = -mask*norm_arr(i,j,k,2)*eb_vel_mag);
                }
             }
           });
       }
     }

     // We make sure to only fill "nghost" ghost cells so we don't accidentally
     // over-write good ghost cell values with unfilled ghost cell values
     IntVect ng_vect(AMREX_D_DECL(nghost,nghost,nghost));
     eb_vel.FillBoundary(0,AMREX_SPACEDIM,ng_vect,gm.periodicity());
}

void
incflo::set_eb_density (int lev, amrex::Real /*time*/, MultiFab& eb_density, int nghost)
{
    Geometry const& gm = Geom(lev);
    eb_density.setVal(0.);

    const auto& factory =
       dynamic_cast<EBFArrayBoxFactory const&>(eb_density.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(eb_density, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       const Box& bx = mfi.tilebox();
       const auto& flagfab      = factory.getMultiEBCellFlagFab()[mfi];

       if (flagfab.getType(bx) == FabType::singlevalued) {
          const auto& flags_arr    = flagfab.const_array();
          const auto& eb_density_arr   = eb_density[mfi].array();
          const auto& norm_arr     = factory.getBndryNormal()[mfi].const_array();

          bool has_normal = m_eb_flow.has_normal;
          GpuArray<Real, AMREX_SPACEDIM> normal{0.};
          if (has_normal) {
             AMREX_D_TERM(
             normal[0] = m_eb_flow.normal[0];,
             normal[1] = m_eb_flow.normal[1];,
             normal[2] = m_eb_flow.normal[2]);
          }
          Real pad = std::numeric_limits<float>::epsilon();
          Real normal_tol = m_eb_flow.normal_tol;
          Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
          Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

          Real eb_flow_density = m_eb_flow.density;

          ParallelFor(bx, [flags_arr,eb_density_arr,norm_arr,has_normal,normal,
                 norm_tol_lo, norm_tol_hi, eb_flow_density]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
             if (flags_arr(i,j,k).isSingleValued()) {
                Real mask = Real(1.0);

                if(has_normal) {
#if (AMREX_SPACEDIM==3)
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1]
                                 + norm_arr(i,j,k,2)*normal[2];
#else
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1];
#endif

                  mask = ((norm_tol_lo <= dotprod) &&
                          (dotprod <= norm_tol_hi)) ? Real(1.0) : Real(0.0);
                }

                eb_density_arr(i,j,k) = mask*eb_flow_density;
             }
           });
       }
     }

     // We make sure to only fill "nghost" ghost cells so we don't accidentally
     // over-write good ghost cell values with unfilled ghost cell values
     IntVect ng_vect(AMREX_D_DECL(nghost,nghost,nghost));
     eb_density.FillBoundary(0,1,ng_vect,gm.periodicity());
}

void
incflo::set_eb_tracer (int lev, amrex::Real /*time*/, MultiFab& eb_tracer, int nghost)
{
    Geometry const& gm = Geom(lev);
    eb_tracer.setVal(0.);

    const auto& factory =
       dynamic_cast<EBFArrayBoxFactory const&>(eb_tracer.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(eb_tracer, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       const Box& bx = mfi.tilebox();
       const auto& flagfab      = factory.getMultiEBCellFlagFab()[mfi];

       if (flagfab.getType(bx) == FabType::singlevalued) {
          const auto& flags_arr    = flagfab.const_array();
          const auto& eb_tracer_arr   = eb_tracer[mfi].array();
          const auto& norm_arr     = factory.getBndryNormal()[mfi].const_array();

          bool has_normal = m_eb_flow.has_normal;
          GpuArray<Real, AMREX_SPACEDIM> normal{0.};
          if (has_normal) {
             AMREX_D_TERM(
             normal[0] = m_eb_flow.normal[0];,
             normal[1] = m_eb_flow.normal[1];,
             normal[2] = m_eb_flow.normal[2]);
          }
          Real pad = std::numeric_limits<float>::epsilon();
          Real normal_tol = m_eb_flow.normal_tol;
          Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
          Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

          int num_trac = m_ntrac;
          // Create a device vector
          Gpu::DeviceVector<Real> eb_flow_tracer_dv;
          for(int n(0); n < num_trac; ++n) {
             eb_flow_tracer_dv.push_back(m_eb_flow.tracer[n]);
          }
          Real* eb_flow_tracer = eb_flow_tracer_dv.data();

          ParallelFor(bx, [flags_arr,eb_tracer_arr,norm_arr,has_normal,normal,
                 norm_tol_lo, norm_tol_hi, num_trac, eb_flow_tracer]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
             if (flags_arr(i,j,k).isSingleValued()) {
                Real mask = Real(1.0);

                if(has_normal) {
#if (AMREX_SPACEDIM==3)
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1]
                                 + norm_arr(i,j,k,2)*normal[2];
#else
                  Real dotprod = norm_arr(i,j,k,0)*normal[0]
                                 + norm_arr(i,j,k,1)*normal[1];
#endif

                  mask = ((norm_tol_lo <= dotprod) &&
                          (dotprod <= norm_tol_hi)) ? Real(1.0) : Real(0.0);
                }

                for(int n(0); n<num_trac; n++) {
                  eb_tracer_arr(i,j,k,n) = mask*eb_flow_tracer[n];
                }
             }
           });
       }
     }

     // We make sure to only fill "nghost" ghost cells so we don't accidentally
     // over-write good ghost cell values with unfilled ghost cell values
     IntVect ng_vect(AMREX_D_DECL(nghost,nghost,nghost));
     eb_tracer.FillBoundary(0,m_ntrac,ng_vect,gm.periodicity());
}
#endif
