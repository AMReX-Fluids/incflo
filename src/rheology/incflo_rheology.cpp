#include <incflo.H>
#include <incflo_derive_K.H>

using namespace amrex;

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real expterm (amrex::Real nu) noexcept
{
    return (nu < Real(1.e-9)) ? (Real(1.0)-Real(0.5)*nu+nu*nu*Real(1.0/6.0)-(nu*nu*nu)*Real(1./24.))
                        : -std::expm1(-nu)/nu;
}

struct NonNewtonianViscosity
{
    incflo::FluidModel fluid_model;
    amrex::Real mu, n_flow, tau_0, eta_0, papa_reg;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real sr) const noexcept {
        switch (fluid_model)
        {
        case incflo::FluidModel::powerlaw:
        {
            return mu * std::pow(sr,n_flow-Real(1.0));
        }
        case incflo::FluidModel::Bingham:
        {
            return mu + tau_0 * expterm(sr/papa_reg) / papa_reg;
        }
        case incflo::FluidModel::HerschelBulkley:
        {
            return (mu*std::pow(sr,n_flow)+tau_0)*expterm(sr/papa_reg)/papa_reg;
        }
        case incflo::FluidModel::deSouzaMendesDutra:
        {
            return (mu*std::pow(sr,n_flow)+tau_0)*expterm(sr*(eta_0/tau_0))*(eta_0/tau_0);
        }
        default:
        {
            return mu;
        }
        };
    }
};

}

void incflo::compute_viscosity (Vector<MultiFab*> const& vel_eta,
                                Vector<MultiFab*> const& rho,
                                Vector<MultiFab*> const& vel,
                                Real time, int nghost)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        compute_viscosity_at_level(lev, vel_eta[lev], rho[lev], vel[lev], geom[lev], time, nghost);
    }
}

#ifdef AMREX_USE_EB
void incflo::compute_viscosity_at_level (int lev,
#else
void incflo::compute_viscosity_at_level (int /*lev*/,
#endif
                                         MultiFab* vel_eta,
                                         MultiFab* /*rho*/,
                                         MultiFab* vel,
                                         Geometry& lev_geom,
                                         Real /*time*/, int nghost)
{
    if (m_fluid_model == FluidModel::Newtonian)
    {
        vel_eta->setVal(m_mu, 0, 1, nghost);
    }
    else
    {
        NonNewtonianViscosity non_newtonian_viscosity;
        non_newtonian_viscosity.fluid_model = m_fluid_model;
        non_newtonian_viscosity.mu = m_mu;
        non_newtonian_viscosity.n_flow = m_n_0;
        non_newtonian_viscosity.tau_0 = m_tau_0;
        non_newtonian_viscosity.eta_0 = m_eta_0;
        non_newtonian_viscosity.papa_reg = m_papa_reg;

#ifdef AMREX_USE_EB
        auto const& fact = EBFactory(lev);
        auto const& flags = fact.getMultiEBCellFlagFab();
#endif

        Real idx = Real(1.0) / lev_geom.CellSize(0);
        Real idy = Real(1.0) / lev_geom.CellSize(1);
#if (AMREX_SPACEDIM == 3)
        Real idz = Real(1.0) / lev_geom.CellSize(2);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vel_eta,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
                Box const& bx = mfi.growntilebox(nghost);
                Array4<Real> const& eta_arr = vel_eta->array(mfi);
                Array4<Real const> const& vel_arr = vel->const_array(mfi);
#ifdef AMREX_USE_EB
                auto const& flag_fab = flags[mfi];
                auto typ = flag_fab.getType(bx);
                if (typ == FabType::covered)
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        eta_arr(i,j,k) = Real(0.0);
                    });
                }
                else if (typ == FabType::singlevalued)
                {
                    auto const& flag_arr = flag_fab.const_array();
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        Real sr = incflo_strainrate_eb(i,j,k,AMREX_D_DECL(idx,idy,idz),vel_arr,flag_arr(i,j,k));
                        eta_arr(i,j,k) = non_newtonian_viscosity(sr);
                    });
                }
                else
#endif
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        Real sr = incflo_strainrate(i,j,k,AMREX_D_DECL(idx,idy,idz),vel_arr);
                        eta_arr(i,j,k) = non_newtonian_viscosity(sr);
                    });
                }
        }
    }
}

void incflo::compute_tracer_diff_coeff (Vector<MultiFab*> const& tra_eta, int nghost)
{
    for (auto mf : tra_eta) {
        for (int n = 0; n < m_ntrac; ++n) {
            mf->setVal(m_mu_s[n], n, 1, nghost);

            // Disk
            if (2976 == m_probtype) {
                // TODO: find a better way to set diff_coeff
                // for (int lev = 0; lev <= finest_level; ++lev) {
                //     auto& ld = *m_leveldata[lev];
                //     Box const& domain = geom[lev].Domain();
                //     auto const& dx = geom[lev].CellSizeArray();
                //     auto const& problo = geom[lev].ProbLoArray();
                //     auto const& probhi = geom[lev].ProbHiArray();
                //     for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
                //     {
                //         const Box& vbx = mfi.validbox();
                //         Array4<Real> const& diff_coeff = mf->array(mfi);
                //         amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                //         {
                //             Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                //             Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                //             Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                //             // Set velocity field inside the plunged objects to 0
                //             Real rz = 3./2;
                //             Real ry = 0.025;
                //             Real disp = m_cur_time*1.;
                //             Real d = disp<(1.+2.*rz)?disp:(1.+2.*rz);
                //             if ((x*x+(z+dx[2]+rz)*(z+dx[2]+rz)<=rz*rz) && (y>=-ry && y<=ry)) {
                //             // if ((x*x+(z-rz+d)*(z-rz+d)<=rz*rz) && (y>=-ry && y<=ry)) {
                //                 diff_coeff(i,j,k) = 0.05547074;
                //             } else diff_coeff(i,j,k) = 0.00016;
                //         });
                //     }
                // }
            }
            // Tcouple
            else if (2977 == m_probtype) {
                // // TODO: find a better way to set diff_coeff
                // for (int lev = 0; lev <= finest_level; ++lev) {
                //     auto& ld = *m_leveldata[lev];
                //     Box const& domain = geom[lev].Domain();
                //     auto const& dx = geom[lev].CellSizeArray();
                //     auto const& problo = geom[lev].ProbLoArray();
                //     auto const& probhi = geom[lev].ProbHiArray();
                //     for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
                //     {
                //         const Box& vbx = mfi.validbox();
                //         Array4<Real> const& diff_coeff = mf->array(mfi);
                //         amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                //         {
                //             Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                //             Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                //             Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                //             // Set velocity field inside the plunged objects to 0
                //             Real r = 0.050;
                //             Real v0 = 0.7408041682457519;
                //             Real acc = -0.04623555291981721;
                //             Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                //             Real d = disp;
                //             if (x*x+y*y+(z-r+dx[2]+d)*(z-r+dx[2]+d)<=r*r) {
                //                 diff_coeff(i,j,k) = 0.00626;
                //             } else diff_coeff(i,j,k) = 0.00025;
                //         });
                //     }
                // }
            }
            // Tcouple with energy
            else if (2978 == m_probtype) {
                // TODO: find a better way to set diff_coeff
                for (int lev = 0; lev <= finest_level; ++lev) {
                    auto& ld = *m_leveldata[lev];
                    Box const& domain = geom[lev].Domain();
                    auto const& dx = geom[lev].CellSizeArray();
                    auto const& problo = geom[lev].ProbLoArray();
                    auto const& probhi = geom[lev].ProbHiArray();
                    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
                    {
                        const Box& vbx = mfi.validbox();
                        Array4<Real> const& diff_coeff = mf->array(mfi);
                        // Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                        // Array4<Real> const& ther_conduct = ld.ther_conduct.array(mfi);
                        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                            Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                            Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                            // Set velocity field inside the plunged objects to 0
                            Real r = 0.050;
                            // // 1 m/s
                            // Real v0 = 0.7408041682457519;
                            // Real acc = -0.04623555291981721;
                            // 500 mm/s
                            Real v0 = 0.485047915352808;
                            Real acc = -0.03131555980546892;
                            Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                            Real d = disp;
                            // Thermal conductivity (converted)
                            Real cv_eth = 1.55;
                            Real cv_tcp = 5.4406;
                            Real k_eth = 0.0003875/cv_eth;
                            Real k_tcp = 0.0340625/cv_tcp;

                            // diff_coeff(i,j,k) = k_tcp/cv_tcp;

                            // Smoothing
                            // Real phi = x*x+y*y+(z+1.5*r+d+dx[2])*(z+1.5*r+d+dx[2]) - r*r;
                            Real phi = x*x+y*y+(z-r+d+dx[2])*(z-r+d+dx[2]) - r*r;
                            Real eps = 2.0*dx[0];
                            if (phi>eps) diff_coeff(i,j,k) = k_eth; // ethane
                            // else diff_coeff(i,j,k) = k_tcp; // tcouple
                            else if (phi<-eps) diff_coeff(i,j,k) = k_tcp; // tcouple
                            else {
                                // 2020 JFM Heaviside function
                                Real xx=phi/eps;
                                diff_coeff(i,j,k) = k_tcp-0.5*(1.+xx+1./M_PI*sin(M_PI*xx))*(k_tcp-k_eth);
                            }
                        });
                    }
                }
            }
            // Disk with energy
            else if (2979 == m_probtype) {
                // TODO: find a better way to set diff_coeff
                for (int lev = 0; lev <= finest_level; ++lev) {
                    auto& ld = *m_leveldata[lev];
                    Box const& domain = geom[lev].Domain();
                    auto const& dx = geom[lev].CellSizeArray();
                    auto const& problo = geom[lev].ProbLoArray();
                    auto const& probhi = geom[lev].ProbHiArray();
                    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
                    {
                        const Box& vbx = mfi.validbox();
                        Array4<Real> const& diff_coeff = mf->array(mfi);
                        // Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                        // Array4<Real> const& ther_conduct = ld.ther_conduct.array(mfi);
                        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                            Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                            Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                            // Set velocity field inside the plunged objects to 0
                            Real rz = 3/2.;
                            Real ry = 0.025;
                            // // 1 m/s
                            // Real v0 = 0.7408041682457519;
                            // Real acc = -0.04623555291981721;
                            // 500 mm/s
                            Real v0 = 0.485047915352808;
                            Real acc = -0.03131555980546892;
                            Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                            Real d = disp;

                            // Thermal conductivity (converted)
                            Real cv_eth = 1.55;
                            Real cv_tcp = 3.76953125; // plunger
                            Real k_eth = 0.0003875/cv_eth;
                            Real k_tcp = 0.000539062/cv_tcp; // plunger

                            if ((x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d)<rz*rz) && (y>-ry && y<ry)) {
                                // tcouple
                                diff_coeff(i,j,k) = k_tcp;
                            } else {
                                // ethane
                                diff_coeff(i,j,k) = k_eth;
                            }
                        });
                    }
                }
            }
            // Grid with energy
            else if (2980 == m_probtype) {
                // TODO: find a better way to set diff_coeff
                for (int lev = 0; lev <= finest_level; ++lev) {
                    auto& ld = *m_leveldata[lev];
                    Box const& domain = geom[lev].Domain();
                    auto const& dx = geom[lev].CellSizeArray();
                    auto const& problo = geom[lev].ProbLoArray();
                    auto const& probhi = geom[lev].ProbHiArray();
                    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
                    {
                        const Box& vbx = mfi.validbox();
                        Array4<Real> const& diff_coeff = mf->array(mfi);
                        // Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                        // Array4<Real> const& ther_conduct = ld.ther_conduct.array(mfi);
                        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                            Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                            Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                            // Set velocity field inside the plunged objects to 0
                            Real rz = 3.05/2.;
                            Real rd = 0.24/2.;
                            Real ry = 0.025/2.;
                            // // 1 m/s
                            // Real v0 = 0.7408041682457519;
                            // Real acc = -0.04623555291981721;
                            // 500 mm/s
                            Real v0 = 0.485047915352808;
                            Real acc = -0.03131555980546892;
                            Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                            Real d = disp;

                            // Thermal conductivity (converted)
                            Real cv_eth = 1.55;
                            Real cv_tcp = 3.76953125; // plunger
                            Real cv_spl = 6.5421875; // sample
                            Real k_eth = 0.0003875/cv_eth;
                            Real k_tcp = 0.000539062/cv_tcp; // plunger
                            Real k_spl = 0.0009344219/cv_spl; // sample

                            Real outer_grid = x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d);
                            Real inner_grid = x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d);
                            Real grid_gap = 0.06;
                            Real zz=z-rz+2*dx[2]+d;

                            if (outer_grid<rz*rz && y*y<ry*ry) {
                                if (inner_grid>(rz-rd)*(rz-rd) ||
                                    (x<ry && x>-ry) ||
                                    (x<ry-1*grid_gap && x>-ry-1*grid_gap) ||
                                    (x<ry-2*grid_gap && x>-ry-2*grid_gap) ||
                                    (x<ry-3*grid_gap && x>-ry-3*grid_gap) ||
                                    (x<ry-4*grid_gap && x>-ry-4*grid_gap) ||
                                    (x<ry-5*grid_gap && x>-ry-5*grid_gap) ||
                                    (x<ry-6*grid_gap && x>-ry-6*grid_gap) ||
                                    (x<ry-7*grid_gap && x>-ry-7*grid_gap) ||
                                    (x<ry-8*grid_gap && x>-ry-8*grid_gap) ||
                                    (x<ry-9*grid_gap && x>-ry-9*grid_gap) ||
                                    (x<ry-10*grid_gap && x>-ry-10*grid_gap) ||
                                    (x<ry-11*grid_gap && x>-ry-11*grid_gap) ||
                                    (x<ry-12*grid_gap && x>-ry-12*grid_gap) ||
                                    (x<ry-13*grid_gap && x>-ry-13*grid_gap) ||
                                    (x<ry-14*grid_gap && x>-ry-14*grid_gap) ||
                                    (x<ry-15*grid_gap && x>-ry-15*grid_gap) ||
                                    (x<ry-16*grid_gap && x>-ry-16*grid_gap) ||
                                    (x<ry-17*grid_gap && x>-ry-17*grid_gap) ||
                                    (x<ry-18*grid_gap && x>-ry-18*grid_gap) ||
                                    (x<ry-19*grid_gap && x>-ry-19*grid_gap) ||
                                    (x<ry-20*grid_gap && x>-ry-20*grid_gap) ||
                                    (x<ry-21*grid_gap && x>-ry-21*grid_gap) ||
                                    (x<ry-22*grid_gap && x>-ry-22*grid_gap) ||
                                    (x<ry-23*grid_gap && x>-ry-23*grid_gap) ||
                                    (x<ry+1*grid_gap && x>-ry+1*grid_gap) ||
                                    (x<ry+2*grid_gap && x>-ry+2*grid_gap) ||
                                    (x<ry+3*grid_gap && x>-ry+3*grid_gap) ||
                                    (x<ry+4*grid_gap && x>-ry+4*grid_gap) ||
                                    (x<ry+5*grid_gap && x>-ry+5*grid_gap) ||
                                    (x<ry+6*grid_gap && x>-ry+6*grid_gap) ||
                                    (x<ry+7*grid_gap && x>-ry+7*grid_gap) ||
                                    (x<ry+8*grid_gap && x>-ry+8*grid_gap) ||
                                    (x<ry+9*grid_gap && x>-ry+9*grid_gap) ||
                                    (x<ry+10*grid_gap && x>-ry+10*grid_gap) ||
                                    (x<ry+11*grid_gap && x>-ry+11*grid_gap) ||
                                    (x<ry+12*grid_gap && x>-ry+12*grid_gap) ||
                                    (x<ry+13*grid_gap && x>-ry+13*grid_gap) ||
                                    (x<ry+14*grid_gap && x>-ry+14*grid_gap) ||
                                    (x<ry+15*grid_gap && x>-ry+15*grid_gap) ||
                                    (x<ry+16*grid_gap && x>-ry+16*grid_gap) ||
                                    (x<ry+17*grid_gap && x>-ry+17*grid_gap) ||
                                    (x<ry+18*grid_gap && x>-ry+18*grid_gap) ||
                                    (x<ry+19*grid_gap && x>-ry+19*grid_gap) ||
                                    (x<ry+20*grid_gap && x>-ry+20*grid_gap) ||
                                    (x<ry+21*grid_gap && x>-ry+21*grid_gap) ||
                                    (x<ry+22*grid_gap && x>-ry+22*grid_gap) ||
                                    (x<ry+23*grid_gap && x>-ry+23*grid_gap) ||
                                    (zz<ry && zz>-ry) ||
                                    (zz<ry-1*grid_gap && zz>-ry-1*grid_gap) ||
                                    (zz<ry-2*grid_gap && zz>-ry-2*grid_gap) ||
                                    (zz<ry-3*grid_gap && zz>-ry-3*grid_gap) ||
                                    (zz<ry-4*grid_gap && zz>-ry-4*grid_gap) ||
                                    (zz<ry-5*grid_gap && zz>-ry-5*grid_gap) ||
                                    (zz<ry-6*grid_gap && zz>-ry-6*grid_gap) ||
                                    (zz<ry-7*grid_gap && zz>-ry-7*grid_gap) ||
                                    (zz<ry-8*grid_gap && zz>-ry-8*grid_gap) ||
                                    (zz<ry-9*grid_gap && zz>-ry-9*grid_gap) ||
                                    (zz<ry-10*grid_gap && zz>-ry-10*grid_gap) ||
                                    (zz<ry-11*grid_gap && zz>-ry-11*grid_gap) ||
                                    (zz<ry-12*grid_gap && zz>-ry-12*grid_gap) ||
                                    (zz<ry-13*grid_gap && zz>-ry-13*grid_gap) ||
                                    (zz<ry-14*grid_gap && zz>-ry-14*grid_gap) ||
                                    (zz<ry-15*grid_gap && zz>-ry-15*grid_gap) ||
                                    (zz<ry-16*grid_gap && zz>-ry-16*grid_gap) ||
                                    (zz<ry-17*grid_gap && zz>-ry-17*grid_gap) ||
                                    (zz<ry-18*grid_gap && zz>-ry-18*grid_gap) ||
                                    (zz<ry-19*grid_gap && zz>-ry-19*grid_gap) ||
                                    (zz<ry-20*grid_gap && zz>-ry-20*grid_gap) ||
                                    (zz<ry-21*grid_gap && zz>-ry-21*grid_gap) ||
                                    (zz<ry-22*grid_gap && zz>-ry-22*grid_gap) ||
                                    (zz<ry-23*grid_gap && zz>-ry-23*grid_gap) ||
                                    (zz<ry+1*grid_gap && zz>-ry+1*grid_gap) ||
                                    (zz<ry+2*grid_gap && zz>-ry+2*grid_gap) ||
                                    (zz<ry+3*grid_gap && zz>-ry+3*grid_gap) ||
                                    (zz<ry+4*grid_gap && zz>-ry+4*grid_gap) ||
                                    (zz<ry+5*grid_gap && zz>-ry+5*grid_gap) ||
                                    (zz<ry+6*grid_gap && zz>-ry+6*grid_gap) ||
                                    (zz<ry+7*grid_gap && zz>-ry+7*grid_gap) ||
                                    (zz<ry+8*grid_gap && zz>-ry+8*grid_gap) ||
                                    (zz<ry+9*grid_gap && zz>-ry+9*grid_gap) ||
                                    (zz<ry+10*grid_gap && zz>-ry+10*grid_gap) ||
                                    (zz<ry+11*grid_gap && zz>-ry+11*grid_gap) ||
                                    (zz<ry+12*grid_gap && zz>-ry+12*grid_gap) ||
                                    (zz<ry+13*grid_gap && zz>-ry+13*grid_gap) ||
                                    (zz<ry+14*grid_gap && zz>-ry+14*grid_gap) ||
                                    (zz<ry+15*grid_gap && zz>-ry+15*grid_gap) ||
                                    (zz<ry+16*grid_gap && zz>-ry+16*grid_gap) ||
                                    (zz<ry+17*grid_gap && zz>-ry+17*grid_gap) ||
                                    (zz<ry+18*grid_gap && zz>-ry+18*grid_gap) ||
                                    (zz<ry+19*grid_gap && zz>-ry+19*grid_gap) ||
                                    (zz<ry+20*grid_gap && zz>-ry+20*grid_gap) ||
                                    (zz<ry+21*grid_gap && zz>-ry+21*grid_gap) ||
                                    (zz<ry+22*grid_gap && zz>-ry+22*grid_gap) ||
                                    (zz<ry+23*grid_gap && zz>-ry+23*grid_gap) ) {
                                    // tcouple
                                    diff_coeff(i,j,k) = k_tcp;
                                }
                            } else {
                                // ethane
                                diff_coeff(i,j,k) = k_eth;
                            }
                            if (inner_grid<0.25*0.25 && (y>ry && y<ry+0.025)) {
                            // if (inner_grid<pow(0.04/2,2) && (y>ry && y<ry+0.01)) {
                                // sample
                                diff_coeff(i,j,k) = k_spl;
                            }
                        });
                    }
                }
            }
        }
    }
}

void incflo::compute_energy_diff_coeff (Vector<MultiFab*> const& energy_eta, int nghost)
{
    for (auto mf : energy_eta) {
        for (int n = 0; n < m_ntrac; ++n) {
            mf->setVal(m_mu_s[n], n, 1, nghost);

            // Disk
            if (2976 == m_probtype) {
                // TODO: find a better way to set diff_coeff
                // for (int lev = 0; lev <= finest_level; ++lev) {
                //     auto& ld = *m_leveldata[lev];
                //     Box const& domain = geom[lev].Domain();
                //     auto const& dx = geom[lev].CellSizeArray();
                //     auto const& problo = geom[lev].ProbLoArray();
                //     auto const& probhi = geom[lev].ProbHiArray();
                //     for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
                //     {
                //         const Box& vbx = mfi.validbox();
                //         Array4<Real> const& diff_coeff = mf->array(mfi);
                //         amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                //         {
                //             Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                //             Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                //             Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                //             // Set velocity field inside the plunged objects to 0
                //             Real rz = 3./2;
                //             Real ry = 0.025;
                //             Real disp = m_cur_time*1.;
                //             Real d = disp<(1.+2.*rz)?disp:(1.+2.*rz);
                //             if ((x*x+(z+dx[2]+rz)*(z+dx[2]+rz)<=rz*rz) && (y>=-ry && y<=ry)) {
                //             // if ((x*x+(z-rz+d)*(z-rz+d)<=rz*rz) && (y>=-ry && y<=ry)) {
                //                 diff_coeff(i,j,k) = 0.05547074;
                //             } else diff_coeff(i,j,k) = 0.00016;
                //         });
                //     }
                // }
            }

            // Tcouple
            else if (2977 == m_probtype) {
                // // TODO: find a better way to set diff_coeff
                // for (int lev = 0; lev <= finest_level; ++lev) {
                //     auto& ld = *m_leveldata[lev];
                //     Box const& domain = geom[lev].Domain();
                //     auto const& dx = geom[lev].CellSizeArray();
                //     auto const& problo = geom[lev].ProbLoArray();
                //     auto const& probhi = geom[lev].ProbHiArray();
                //     for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
                //     {
                //         const Box& vbx = mfi.validbox();
                //         Array4<Real> const& diff_coeff = mf->array(mfi);
                //         amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                //         {
                //             Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                //             Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                //             Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                //             // Set velocity field inside the plunged objects to 0
                //             Real r = 0.050;
                //             Real v0 = 0.7408041682457519;
                //             Real acc = -0.04623555291981721;
                //             Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                //             Real d = disp;
                //             if (x*x+y*y+(z-r+dx[2]+d)*(z-r+dx[2]+d)<=r*r) {
                //                 diff_coeff(i,j,k) = 0.00626;
                //             } else diff_coeff(i,j,k) = 0.00025;
                //         });
                //     }
                // }
            }

            // Tcouple with energy
            else if (2978 == m_probtype) {
                // TODO: find a better way to set diff_coeff
                for (int lev = 0; lev <= finest_level; ++lev) {
                    auto& ld = *m_leveldata[lev];
                    Box const& domain = geom[lev].Domain();
                    auto const& dx = geom[lev].CellSizeArray();
                    auto const& problo = geom[lev].ProbLoArray();
                    auto const& probhi = geom[lev].ProbHiArray();
                    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
                    {
                        const Box& vbx = mfi.validbox();
                        Array4<Real> const& diff_coeff = mf->array(mfi);
                        // Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                        // Array4<Real> const& ther_conduct = ld.ther_conduct.array(mfi);
                        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                            Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                            Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                            // Set velocity field inside the plunged objects to 0
                            Real r = 0.050;
                            // // 1 m/s
                            // Real v0 = 0.7408041682457519;
                            // Real acc = -0.04623555291981721;
                            // 500 mm/s
                            Real v0 = 0.485047915352808;
                            Real acc = -0.03131555980546892;
                            Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                            Real d = disp;
                            // Thermal conductivity (converted)
                            Real cv_eth = 1.55;
                            Real cv_tcp = 5.4406;
                            Real k_eth = 0.0003875;//cv_eth;
                            Real k_tcp = 0.0340625;//cv_tcp;

                            // diff_coeff(i,j,k) = k_tcp/cv_tcp;

                            // Smoothing
                            // Real phi = x*x+y*y+(z+1.5*r+d+dx[2])*(z+1.5*r+d+dx[2]) - r*r;
                            Real phi = x*x+y*y+(z-r+d+dx[2])*(z-r+d+dx[2]) - r*r;
                            Real eps = 2.0*dx[0];
                            if (phi>eps) diff_coeff(i,j,k) = k_eth; // ethane
                            // else diff_coeff(i,j,k) = k_tcp; // tcouple
                            else if (phi<-eps) diff_coeff(i,j,k) = k_tcp; // tcouple
                            else {
                                // 2020 JFM Heaviside function
                                Real xx=phi/eps;
                                diff_coeff(i,j,k) = k_tcp-0.5*(1.+xx+1./M_PI*sin(M_PI*xx))*(k_tcp-k_eth);
                            }
                        });
                    }
                }
            }

            // Disk with energy
            else if (2979 == m_probtype) {
                // TODO: find a better way to set diff_coeff
                for (int lev = 0; lev <= finest_level; ++lev) {
                    auto& ld = *m_leveldata[lev];
                    Box const& domain = geom[lev].Domain();
                    auto const& dx = geom[lev].CellSizeArray();
                    auto const& problo = geom[lev].ProbLoArray();
                    auto const& probhi = geom[lev].ProbHiArray();
                    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
                    {
                        const Box& vbx = mfi.validbox();
                        Array4<Real> const& diff_coeff = mf->array(mfi);
                        // Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                        // Array4<Real> const& ther_conduct = ld.ther_conduct.array(mfi);
                        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                            Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                            Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                            // Set velocity field inside the plunged objects to 0
                            Real rz = 3/2.;
                            Real ry = 0.025;
                            // // 1 m/s
                            // Real v0 = 0.7408041682457519;
                            // Real acc = -0.04623555291981721;
                            // 500 mm/s
                            Real v0 = 0.485047915352808;
                            Real acc = -0.03131555980546892;
                            Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                            Real d = disp;

                            // Thermal conductivity (converted)
                            Real cv_eth = 1.55;
                            Real cv_tcp = 3.76953125; //plunger
                            Real k_eth = 0.0003875;
                            Real k_tcp = 0.0005390625; // plunger

                            if ((x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d)<rz*rz) && (y>-ry && y<ry)) {
                                // tcouple
                                diff_coeff(i,j,k) = k_tcp;
                            } else {
                                // ethane
                                diff_coeff(i,j,k) = k_eth;
                            }
                        });
                    }
                }
            }

            // Grid with energy
            else if (2980 == m_probtype) {
                // TODO: find a better way to set diff_coeff
                for (int lev = 0; lev <= finest_level; ++lev) {
                    auto& ld = *m_leveldata[lev];
                    Box const& domain = geom[lev].Domain();
                    auto const& dx = geom[lev].CellSizeArray();
                    auto const& problo = geom[lev].ProbLoArray();
                    auto const& probhi = geom[lev].ProbHiArray();
                    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
                    {
                        const Box& vbx = mfi.validbox();
                        Array4<Real> const& diff_coeff = mf->array(mfi);
                        // Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                        // Array4<Real> const& ther_conduct = ld.ther_conduct.array(mfi);
                        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                            Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                            Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                            // Set velocity field inside the plunged objects to 0
                            Real rz = 3.05/2.;
                            Real rd = 0.24/2.;
                            Real ry = 0.025/2.;
                            // // 1 m/s
                            // Real v0 = 0.7408041682457519;
                            // Real acc = -0.04623555291981721;
                            // 500 mm/s
                            Real v0 = 0.485047915352808;
                            Real acc = -0.03131555980546892;
                            Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                            Real d = disp;

                            // Thermal conductivity (converted)
                            Real cv_eth = 1.55;
                            Real cv_tcp = 3.76953125; //plunger
                            Real cv_spl = 6.5421875; // sample
                            Real k_eth = 0.0003875;
                            Real k_tcp = 0.0005390625; // plunger
                            Real k_spl = 0.0009344219; // sample

                            Real outer_grid = x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d);
                            Real inner_grid = x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d);
                            Real grid_gap = 0.06;
                            Real zz=z-rz+2*dx[2]+d;

                            if (outer_grid<rz*rz && y*y<ry*ry) {
                                if (inner_grid>(rz-rd)*(rz-rd) ||
                                    (x<ry && x>-ry) ||
                                    (x<ry-1*grid_gap && x>-ry-1*grid_gap) ||
                                    (x<ry-2*grid_gap && x>-ry-2*grid_gap) ||
                                    (x<ry-3*grid_gap && x>-ry-3*grid_gap) ||
                                    (x<ry-4*grid_gap && x>-ry-4*grid_gap) ||
                                    (x<ry-5*grid_gap && x>-ry-5*grid_gap) ||
                                    (x<ry-6*grid_gap && x>-ry-6*grid_gap) ||
                                    (x<ry-7*grid_gap && x>-ry-7*grid_gap) ||
                                    (x<ry-8*grid_gap && x>-ry-8*grid_gap) ||
                                    (x<ry-9*grid_gap && x>-ry-9*grid_gap) ||
                                    (x<ry-10*grid_gap && x>-ry-10*grid_gap) ||
                                    (x<ry-11*grid_gap && x>-ry-11*grid_gap) ||
                                    (x<ry-12*grid_gap && x>-ry-12*grid_gap) ||
                                    (x<ry-13*grid_gap && x>-ry-13*grid_gap) ||
                                    (x<ry-14*grid_gap && x>-ry-14*grid_gap) ||
                                    (x<ry-15*grid_gap && x>-ry-15*grid_gap) ||
                                    (x<ry-16*grid_gap && x>-ry-16*grid_gap) ||
                                    (x<ry-17*grid_gap && x>-ry-17*grid_gap) ||
                                    (x<ry-18*grid_gap && x>-ry-18*grid_gap) ||
                                    (x<ry-19*grid_gap && x>-ry-19*grid_gap) ||
                                    (x<ry-20*grid_gap && x>-ry-20*grid_gap) ||
                                    (x<ry-21*grid_gap && x>-ry-21*grid_gap) ||
                                    (x<ry-22*grid_gap && x>-ry-22*grid_gap) ||
                                    (x<ry-23*grid_gap && x>-ry-23*grid_gap) ||
                                    (x<ry+1*grid_gap && x>-ry+1*grid_gap) ||
                                    (x<ry+2*grid_gap && x>-ry+2*grid_gap) ||
                                    (x<ry+3*grid_gap && x>-ry+3*grid_gap) ||
                                    (x<ry+4*grid_gap && x>-ry+4*grid_gap) ||
                                    (x<ry+5*grid_gap && x>-ry+5*grid_gap) ||
                                    (x<ry+6*grid_gap && x>-ry+6*grid_gap) ||
                                    (x<ry+7*grid_gap && x>-ry+7*grid_gap) ||
                                    (x<ry+8*grid_gap && x>-ry+8*grid_gap) ||
                                    (x<ry+9*grid_gap && x>-ry+9*grid_gap) ||
                                    (x<ry+10*grid_gap && x>-ry+10*grid_gap) ||
                                    (x<ry+11*grid_gap && x>-ry+11*grid_gap) ||
                                    (x<ry+12*grid_gap && x>-ry+12*grid_gap) ||
                                    (x<ry+13*grid_gap && x>-ry+13*grid_gap) ||
                                    (x<ry+14*grid_gap && x>-ry+14*grid_gap) ||
                                    (x<ry+15*grid_gap && x>-ry+15*grid_gap) ||
                                    (x<ry+16*grid_gap && x>-ry+16*grid_gap) ||
                                    (x<ry+17*grid_gap && x>-ry+17*grid_gap) ||
                                    (x<ry+18*grid_gap && x>-ry+18*grid_gap) ||
                                    (x<ry+19*grid_gap && x>-ry+19*grid_gap) ||
                                    (x<ry+20*grid_gap && x>-ry+20*grid_gap) ||
                                    (x<ry+21*grid_gap && x>-ry+21*grid_gap) ||
                                    (x<ry+22*grid_gap && x>-ry+22*grid_gap) ||
                                    (x<ry+23*grid_gap && x>-ry+23*grid_gap) ||
                                    (zz<ry && zz>-ry) ||
                                    (zz<ry-1*grid_gap && zz>-ry-1*grid_gap) ||
                                    (zz<ry-2*grid_gap && zz>-ry-2*grid_gap) ||
                                    (zz<ry-3*grid_gap && zz>-ry-3*grid_gap) ||
                                    (zz<ry-4*grid_gap && zz>-ry-4*grid_gap) ||
                                    (zz<ry-5*grid_gap && zz>-ry-5*grid_gap) ||
                                    (zz<ry-6*grid_gap && zz>-ry-6*grid_gap) ||
                                    (zz<ry-7*grid_gap && zz>-ry-7*grid_gap) ||
                                    (zz<ry-8*grid_gap && zz>-ry-8*grid_gap) ||
                                    (zz<ry-9*grid_gap && zz>-ry-9*grid_gap) ||
                                    (zz<ry-10*grid_gap && zz>-ry-10*grid_gap) ||
                                    (zz<ry-11*grid_gap && zz>-ry-11*grid_gap) ||
                                    (zz<ry-12*grid_gap && zz>-ry-12*grid_gap) ||
                                    (zz<ry-13*grid_gap && zz>-ry-13*grid_gap) ||
                                    (zz<ry-14*grid_gap && zz>-ry-14*grid_gap) ||
                                    (zz<ry-15*grid_gap && zz>-ry-15*grid_gap) ||
                                    (zz<ry-16*grid_gap && zz>-ry-16*grid_gap) ||
                                    (zz<ry-17*grid_gap && zz>-ry-17*grid_gap) ||
                                    (zz<ry-18*grid_gap && zz>-ry-18*grid_gap) ||
                                    (zz<ry-19*grid_gap && zz>-ry-19*grid_gap) ||
                                    (zz<ry-20*grid_gap && zz>-ry-20*grid_gap) ||
                                    (zz<ry-21*grid_gap && zz>-ry-21*grid_gap) ||
                                    (zz<ry-22*grid_gap && zz>-ry-22*grid_gap) ||
                                    (zz<ry-23*grid_gap && zz>-ry-23*grid_gap) ||
                                    (zz<ry+1*grid_gap && zz>-ry+1*grid_gap) ||
                                    (zz<ry+2*grid_gap && zz>-ry+2*grid_gap) ||
                                    (zz<ry+3*grid_gap && zz>-ry+3*grid_gap) ||
                                    (zz<ry+4*grid_gap && zz>-ry+4*grid_gap) ||
                                    (zz<ry+5*grid_gap && zz>-ry+5*grid_gap) ||
                                    (zz<ry+6*grid_gap && zz>-ry+6*grid_gap) ||
                                    (zz<ry+7*grid_gap && zz>-ry+7*grid_gap) ||
                                    (zz<ry+8*grid_gap && zz>-ry+8*grid_gap) ||
                                    (zz<ry+9*grid_gap && zz>-ry+9*grid_gap) ||
                                    (zz<ry+10*grid_gap && zz>-ry+10*grid_gap) ||
                                    (zz<ry+11*grid_gap && zz>-ry+11*grid_gap) ||
                                    (zz<ry+12*grid_gap && zz>-ry+12*grid_gap) ||
                                    (zz<ry+13*grid_gap && zz>-ry+13*grid_gap) ||
                                    (zz<ry+14*grid_gap && zz>-ry+14*grid_gap) ||
                                    (zz<ry+15*grid_gap && zz>-ry+15*grid_gap) ||
                                    (zz<ry+16*grid_gap && zz>-ry+16*grid_gap) ||
                                    (zz<ry+17*grid_gap && zz>-ry+17*grid_gap) ||
                                    (zz<ry+18*grid_gap && zz>-ry+18*grid_gap) ||
                                    (zz<ry+19*grid_gap && zz>-ry+19*grid_gap) ||
                                    (zz<ry+20*grid_gap && zz>-ry+20*grid_gap) ||
                                    (zz<ry+21*grid_gap && zz>-ry+21*grid_gap) ||
                                    (zz<ry+22*grid_gap && zz>-ry+22*grid_gap) ||
                                    (zz<ry+23*grid_gap && zz>-ry+23*grid_gap) ) {
                                    // tcouple
                                    diff_coeff(i,j,k) = k_tcp;
                                }
                            } else {
                                // ethane
                                diff_coeff(i,j,k) = k_eth;
                            }
                            if (inner_grid<0.25*0.25 && (y>ry && y<ry+0.025)) {
                                // sample
                                diff_coeff(i,j,k) = k_spl;
                            }
                        });
                    }
                }
            }
        }
    }
}
