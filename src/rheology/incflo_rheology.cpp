#include <incflo.H>
#include <incflo_derive_K.H>

using namespace amrex;

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real expterm (amrex::Real nu) noexcept
{
    return (nu < 1.e-9) ? (1.0-0.5*nu+nu*nu*(1.0/6.0)-(nu*nu*nu)*(1./24.))
                        : -std::expm1(-nu)/nu;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real get_concentration (amrex::Real dens, amrex::Real rho0, amrex::Real rho1) noexcept
{
    return ((rho1/dens) - 1.0)/((rho1/rho0) - 1.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real mixture_viscosity (amrex::Real conc, amrex::Real visc0, amrex::Real visc1) noexcept
{
    return visc1/(1.0 + conc*((visc1/visc0) - 1.0));
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
            return mu * std::pow(sr,n_flow-1.0);
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

struct ViscosityVOF
{
    amrex::Vector<incflo::FluidModel> fluid_model{{incflo::FluidModel::Newtonian, incflo::FluidModel::Newtonian}};
    amrex::Vector<amrex::Real> rho{{1.0, 1.0}};
    amrex::Vector<amrex::Real> mu{{1.0, 1.0}}; 
    amrex::Vector<amrex::Real> n_flow{{0.0, 0.0}}; 
    amrex::Vector<amrex::Real> tau_0{{0.0, 0.0}}; 
    amrex::Vector<amrex::Real> eta_0{{0.0, 0.0}}; 
    amrex::Vector<amrex::Real> papa_reg{{0.0, 0.0}}; 

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real sr, amrex::Real density) const noexcept {
        amrex::Real visc0;
        if (fluid_model[0] ==  incflo::FluidModel::powerlaw)
        {
            visc0 =  mu[0] * std::pow(sr,n_flow[0]-1.0);
        }
        else if (fluid_model[0] == incflo::FluidModel::Bingham)
        {
            visc0 =  mu[0] + tau_0[0] * expterm(sr/papa_reg[0]) / papa_reg[0];
        }
        else if (fluid_model[0] == incflo::FluidModel::HerschelBulkley)
        {
            visc0 =  (mu[0]*std::pow(sr,n_flow[0])+tau_0[0])*expterm(sr/papa_reg[0])/papa_reg[0];
        }
        else if (fluid_model[0] == incflo::FluidModel::deSouzaMendesDutra)
        {
            visc0 =  (mu[0]*std::pow(sr,n_flow[0])+tau_0[0])*
                      expterm(sr*(eta_0[0]/tau_0[0]))*(eta_0[0]/tau_0[0]);
        }
        else
        {
            visc0 =  mu[0];
        }

        amrex::Real visc1;
        if (fluid_model[1] == incflo::FluidModel::powerlaw)
        {
            visc1 =  mu[1] * std::pow(sr,n_flow[1]-1.0);
        }
        else if (fluid_model[1] == incflo::FluidModel::Bingham)
        {
            visc1 =  mu[1] + tau_0[1] * expterm(sr/papa_reg[1]) / papa_reg[1];
        }
        else if (fluid_model[1] == incflo::FluidModel::HerschelBulkley)
        {
            visc1 =  (mu[1]*std::pow(sr,n_flow[1])+tau_0[1])*expterm(sr/papa_reg[1])/papa_reg[1];
        }
        else if (fluid_model[1] == incflo::FluidModel::deSouzaMendesDutra)
        {
            visc1 =  (mu[1]*std::pow(sr,n_flow[1])+tau_0[1])*
                      expterm(sr*(eta_0[1]/tau_0[1]))*(eta_0[1]/tau_0[1]);
        }
        else
        {
            visc1 =  mu[1];
        }

        amrex::Real conc = get_concentration(density,rho[0],rho[1]);
        return mixture_viscosity(conc,visc0,visc1);

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
                                         MultiFab* rho,
                                         MultiFab* vel,
                                         Geometry& lev_geom,
                                         Real /*time*/, int nghost)
{
    if ((m_fluid_model == FluidModel::Newtonian) and (!m_do_vof))
    {
        vel_eta->setVal(m_mu, 0, 1, nghost);
    }
    else
    {
        ViscosityVOF viscosity_vof;
        NonNewtonianViscosity non_newtonian_viscosity;
        if (m_do_vof) {
            viscosity_vof.fluid_model[0] = m_fluid_vof[0].fluid_model;
            viscosity_vof.rho[0] = m_fluid_vof[0].rho;
            viscosity_vof.mu[0] = m_fluid_vof[0].mu;
            viscosity_vof.n_flow[0] = m_fluid_vof[0].n_0;
            viscosity_vof.tau_0[0] = m_fluid_vof[0].tau_0;
            viscosity_vof.eta_0[0] = m_fluid_vof[0].eta_0;
            viscosity_vof.papa_reg[0] = m_fluid_vof[0].papa_reg;

            viscosity_vof.fluid_model[1] = m_fluid_vof[1].fluid_model;
            viscosity_vof.rho[1] = m_fluid_vof[1].rho;
            viscosity_vof.mu[1] = m_fluid_vof[1].mu;
            viscosity_vof.n_flow[1] = m_fluid_vof[1].n_0;
            viscosity_vof.tau_0[1] = m_fluid_vof[1].tau_0;
            viscosity_vof.eta_0[1] = m_fluid_vof[1].eta_0;
            viscosity_vof.papa_reg[1] = m_fluid_vof[1].papa_reg;
        }
        else {
            non_newtonian_viscosity.fluid_model = m_fluid_model;
            non_newtonian_viscosity.mu = m_mu;
            non_newtonian_viscosity.n_flow = m_n_0;
            non_newtonian_viscosity.tau_0 = m_tau_0;
            non_newtonian_viscosity.eta_0 = m_eta_0;
            non_newtonian_viscosity.papa_reg = m_papa_reg;
        }

#ifdef AMREX_USE_EB
        auto const& fact = EBFactory(lev);
        auto const& flags = fact.getMultiEBCellFlagFab();
#endif

        Real idx = 1.0 / lev_geom.CellSize(0);
        Real idy = 1.0 / lev_geom.CellSize(1);
#if (AMREX_SPACEDIM == 3)
        Real idz = 1.0 / lev_geom.CellSize(2);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vel_eta,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
                Box const& bx = mfi.growntilebox(nghost);
                Array4<Real> const& eta_arr = vel_eta->array(mfi);
                Array4<Real const> const& vel_arr = vel->const_array(mfi);
                Array4<Real const> const& rho_arr = rho->const_array(mfi);
#ifdef AMREX_USE_EB
                auto const& flag_fab = flags[mfi];
                auto typ = flag_fab.getType(bx);
                if (typ == FabType::covered)
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        eta_arr(i,j,k) = 0.0;
                    });
                }
                else if (typ == FabType::singlevalued)
                {
                    auto const& flag_arr = flag_fab.const_array();
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        Real sr = incflo_strainrate_eb(i,j,k,AMREX_D_DECL(idx,idy,idz),vel_arr,flag_arr(i,j,k));
                        Real dens = rho_arr(i,j,k);
                        if (m_do_vof) eta_arr(i,j,k) = viscosity_vof(sr,dens);
                        else eta_arr(i,j,k) = non_newtonian_viscosity(sr);
                    });
                }
                else
#endif
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        Real sr = incflo_strainrate(i,j,k,AMREX_D_DECL(idx,idy,idz),vel_arr);
                        Real dens = rho_arr(i,j,k);
                        if (m_do_vof) eta_arr(i,j,k) = viscosity_vof(sr,dens);
                        else eta_arr(i,j,k) = non_newtonian_viscosity(sr);
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
        }
    }
}
