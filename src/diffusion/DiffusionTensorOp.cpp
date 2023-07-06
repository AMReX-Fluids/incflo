#include <incflo.H>
#include <DiffusionTensorOp.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB_Redistribution.H>
#include <memory>
#endif

using namespace amrex;

DiffusionTensorOp::DiffusionTensorOp (incflo* a_incflo)
    : m_incflo(a_incflo)
{
    define(m_incflo->m_cur_time);
}

DiffusionTensorOp::DiffusionTensorOp (incflo* a_incflo, Real time)
    : m_incflo(a_incflo)
{
    define(time);
}

void DiffusionTensorOp::define (Real time)
{
    readParameters();

    int finest_level = m_incflo->finestLevel();

    LPInfo info_solve;
    info_solve.setMaxCoarseningLevel(m_mg_max_coarsening_level);
    LPInfo info_apply;
    info_apply.setMaxCoarseningLevel(0);
#ifdef AMREX_USE_EB
    if (!m_incflo->EBFactory(0, time).isAllRegular())
    {
        Vector<EBFArrayBoxFactory const*> ebfact;
        for (int lev = 0; lev <= finest_level; ++lev) {
            ebfact.push_back(&(m_incflo->EBFactory(lev,time)));
        }

        if (m_incflo->useTensorSolve())
        {
            m_eb_solve_op = std::make_unique<MLEBTensorOp>(m_incflo->Geom(0,finest_level),
                                                 m_incflo->boxArray(0,finest_level),
                                                 m_incflo->DistributionMap(0,finest_level),
                                                 info_solve, ebfact);
            m_eb_solve_op->setMaxOrder(m_mg_maxorder);
            m_eb_solve_op->setDomainBC(m_incflo->get_diffuse_tensor_bc(Orientation::low),
                                       m_incflo->get_diffuse_tensor_bc(Orientation::high));
        }

        if (m_incflo->need_divtau() || m_incflo->useTensorCorrection())
        {
            m_eb_apply_op = std::make_unique<MLEBTensorOp>(m_incflo->Geom(0,finest_level),
                                                 m_incflo->boxArray(0,finest_level),
                                                 m_incflo->DistributionMap(0,finest_level),
                                                 info_apply, ebfact);
            m_eb_apply_op->setMaxOrder(m_mg_maxorder);
            m_eb_apply_op->setDomainBC(m_incflo->get_diffuse_tensor_bc(Orientation::low),
                                       m_incflo->get_diffuse_tensor_bc(Orientation::high));
        }
    }
    else
#endif
    {
        if (m_incflo->useTensorSolve())
        {
            m_reg_solve_op = std::make_unique<MLTensorOp>(m_incflo->Geom(0,finest_level),
                                                m_incflo->boxArray(0,finest_level),
                                                m_incflo->DistributionMap(0,finest_level),
                                                info_solve);
            m_reg_solve_op->setMaxOrder(m_mg_maxorder);
            m_reg_solve_op->setDomainBC(m_incflo->get_diffuse_tensor_bc(Orientation::low),
                                        m_incflo->get_diffuse_tensor_bc(Orientation::high));
        }

        if (m_incflo->need_divtau() || m_incflo->useTensorCorrection())
        {
            m_reg_apply_op = std::make_unique<MLTensorOp>(m_incflo->Geom(0,finest_level),
                                                m_incflo->boxArray(0,finest_level),
                                                m_incflo->DistributionMap(0,finest_level),
                                                info_apply);
            m_reg_apply_op->setMaxOrder(m_mg_maxorder);
            m_reg_apply_op->setDomainBC(m_incflo->get_diffuse_tensor_bc(Orientation::low),
                                        m_incflo->get_diffuse_tensor_bc(Orientation::high));
        }
    }
}

void
DiffusionTensorOp::readParameters ()
{
    ParmParse pp("tensor_diffusion");

    pp.query("verbose", m_verbose);
    pp.query("mg_verbose", m_mg_verbose);
    pp.query("mg_bottom_verbose", m_mg_bottom_verbose);
    pp.query("mg_max_iter", m_mg_max_iter);
    pp.query("mg_bottom_maxiter", m_mg_bottom_maxiter);
    pp.query("mg_max_fmg_iter", m_mg_max_fmg_iter);
    pp.query("mg_max_coarsening_level", m_mg_max_coarsening_level);
    pp.query("mg_maxorder", m_mg_maxorder);
    pp.query("mg_rtol", m_mg_rtol);
    pp.query("mg_atol", m_mg_atol);
    pp.query("bottom_solver", m_bottom_solver);

    pp.query("num_pre_smooth", m_num_pre_smooth);
    pp.query("num_post_smooth", m_num_post_smooth);
}

void
DiffusionTensorOp::diffuse_velocity (Vector<MultiFab*> const& velocity,
                                     Vector<MultiFab*> const& density,
                                     Vector<MultiFab const*> const& eta,
                                     Real dt)
{
    //
    //      alpha a - beta div ( b grad )   <--->   rho - dt div ( mu grad )
    //
    // So the constants and variable coefficients are:
    //
    //      alpha: 1
    //      beta: dt
    //      a: rho
    //      b: mu

    if (m_verbose > 0) {
        amrex::Print() << "Diffusing velocity components all together..." << std::endl;
    }

    const int finest_level = m_incflo->finestLevel();

#ifdef AMREX_USE_EB
    if (m_eb_solve_op)
    {
        // For when we use the stencil for centroid values
        // m_eb_solve_op->setPhiOnCentroid();

        m_eb_solve_op->setScalars(1.0, dt);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_eb_solve_op->setACoeffs(lev, *density[lev]);

            Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_velocity_eta_to_faces(lev, *eta[lev]);

            m_eb_solve_op->setShearViscosity(lev, GetArrOfConstPtrs(b), MLMG::Location::FaceCentroid);

            if (m_incflo->hasEBFlow()) {
               m_eb_solve_op->setEBShearViscosityWithInflow(lev, *eta[lev], *(m_incflo->get_velocity_eb()[lev]));
            } else {
               m_eb_solve_op->setEBShearViscosity(lev, *eta[lev]);
            }
        }
    }
    else
#endif
    {
        m_reg_solve_op->setScalars(1.0, dt);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_reg_solve_op->setACoeffs(lev, *density[lev]);
            Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_velocity_eta_to_faces(lev, *eta[lev]);
            m_reg_solve_op->setShearViscosity(lev, GetArrOfConstPtrs(b));
        }
    }

    Vector<MultiFab> rhs(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rhs[lev].define(velocity[lev]->boxArray(),
                        velocity[lev]->DistributionMap(), AMREX_SPACEDIM, 0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& rhs_a = rhs[lev].array(mfi);
            Array4<Real const> const& vel_a = velocity[lev]->const_array(mfi);
            Array4<Real const> const& rho_a = density[lev]->const_array(mfi);
            amrex::ParallelFor(bx,AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                rhs_a(i,j,k,n) = rho_a(i,j,k) * vel_a(i,j,k,n);
            });
        }

#ifdef AMREX_USE_EB
        if (m_eb_solve_op) {
            m_eb_solve_op->setLevelBC(lev, velocity[lev]);
        } else
#endif
        {
            m_reg_solve_op->setLevelBC(lev, velocity[lev]);
        }
    }

#ifdef AMREX_USE_EB
    MLMG mlmg(m_eb_solve_op ? static_cast<MLLinOp&>(*m_eb_solve_op)
              :               static_cast<MLLinOp&>(*m_reg_solve_op));
#else
    MLMG mlmg(*m_reg_solve_op);
#endif

    // The default bottom solver is BiCG
    if (m_bottom_solver == "smoother")
    {
        mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if (m_bottom_solver == "hypre")
    {
        mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
    }
    // Maximum iterations for MultiGrid / ConjugateGradients
    mlmg.setMaxIter(m_mg_max_iter);
    mlmg.setMaxFmgIter(m_mg_max_fmg_iter);
    mlmg.setBottomMaxIter(m_mg_bottom_maxiter);

    // Verbosity for MultiGrid / ConjugateGradients
    mlmg.setVerbose(m_mg_verbose);
    mlmg.setBottomVerbose(m_mg_bottom_verbose);

    mlmg.setPreSmooth(m_num_pre_smooth);
    mlmg.setPostSmooth(m_num_post_smooth);

    mlmg.solve(velocity, GetVecOfConstPtrs(rhs), m_mg_rtol, m_mg_atol);
}

void DiffusionTensorOp::compute_divtau (Vector<MultiFab*> const& a_divtau,
                                        Vector<MultiFab const*> const& a_velocity,
                                        Vector<MultiFab const*> const& a_density,
                                        Vector<MultiFab const*> const& a_eta)
{
    BL_PROFILE("DiffusionTensorOp::compute_divtau");

    int finest_level = m_incflo->finestLevel();

    // FIXME? This copy is needed because MLMG::apply(... , const Vector<MF*>&), so
    // can't pass Vector<MF const*>
    Vector<MultiFab> velocity(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        velocity[lev].define(a_velocity[lev]->boxArray(),
                             a_velocity[lev]->DistributionMap(),
                             AMREX_SPACEDIM, 1, MFInfo(),
                             a_velocity[lev]->Factory());
        MultiFab::Copy(velocity[lev], *a_velocity[lev], 0, 0, AMREX_SPACEDIM, 1);
    }

#ifdef AMREX_USE_EB
    if (m_eb_apply_op)
    {
        // We want to return div (mu grad)) phi
        m_eb_apply_op->setScalars(0.0, -1.0);

        // For when we use the stencil for centroid values
        // m_eb_apply_op->setPhiOnCentroid();

        for (int lev = 0; lev <= finest_level; ++lev) {
            m_eb_apply_op->setACoeffs(lev, *a_density[lev]);

            Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_velocity_eta_to_faces(lev, *a_eta[lev]);

            m_eb_apply_op->setShearViscosity(lev, GetArrOfConstPtrs(b), MLMG::Location::FaceCentroid);

            if (m_incflo->hasEBFlow()) {
               m_eb_apply_op->setEBShearViscosityWithInflow(lev, *a_eta[lev], *(m_incflo->get_velocity_eb()[lev]));
            } else {
               m_eb_apply_op->setEBShearViscosity(lev, *a_eta[lev]);
            }
            m_eb_apply_op->setLevelBC(lev, &velocity[lev]);
        }


#ifdef AMREX_USE_MOVING_EB
        //
        // For moving EB, don't redistribute yet.
        //
        MLMG mlmg(*m_eb_apply_op);
        mlmg.apply(a_divtau, GetVecOfPtrs(velocity));

#else
        //
        // Need a temporary to redistribute this term before potential use in
        // an implicit solve
        //
        Vector<MultiFab> divtau_tmp(finest_level+1);
        int tmp_comp = (m_incflo->m_redistribution_type == "StateRedist") ? 3 : 2;
        for (int lev = 0; lev <= finest_level; ++lev) {
            divtau_tmp[lev].define(a_divtau[lev]->boxArray(),
                                   a_divtau[lev]->DistributionMap(),
                                   AMREX_SPACEDIM, tmp_comp, MFInfo(),
                                   a_divtau[lev]->Factory());
            divtau_tmp[lev].setVal(0.0);
        }

        MLMG mlmg(*m_eb_apply_op);
        mlmg.apply(GetVecOfPtrs(divtau_tmp), GetVecOfPtrs(velocity));

        for(int lev = 0; lev <= finest_level; lev++)
        {
            // Flux redistribution
            amrex::single_level_redistribute( divtau_tmp[lev], *a_divtau[lev], 0, AMREX_SPACEDIM, m_incflo->Geom(lev));
            //
            // If we want to allow option of SRD, use incflo::redistribute_term.
            //
            // auto const& bc = m_incflo->get_velocity_bcrec_device_ptr();
            // m_incflo->redistribute_term(*a_divtau[lev], divtau_tmp[lev], *a_velocity[lev],
            //                          bc, lev, Array4<Real const>{});
        }
#endif

    }
    else
#endif
    {
        // We want to return div (mu grad)) phi
        m_reg_apply_op->setScalars(0.0, -1.0);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_reg_apply_op->setACoeffs(lev, *a_density[lev]);
            Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_velocity_eta_to_faces(lev, *a_eta[lev]);
            m_reg_apply_op->setShearViscosity(lev, GetArrOfConstPtrs(b));
            m_reg_apply_op->setLevelBC(lev, &velocity[lev]);
        }

        MLMG mlmg(*m_reg_apply_op);
        mlmg.apply(a_divtau, GetVecOfPtrs(velocity));
    }

    bool advect_momentum = m_incflo->AdvectMomentum();
    if (!advect_momentum) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (MFIter mfi(*a_divtau[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& divtau_arr = a_divtau[lev]->array(mfi);
                Array4<Real const> const& rho_arr = a_density[lev]->const_array(mfi);
                    amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhoinv = Real(1.0)/rho_arr(i,j,k);
                    AMREX_D_TERM(divtau_arr(i,j,k,0) *= rhoinv;,
                                 divtau_arr(i,j,k,1) *= rhoinv;,
                                 divtau_arr(i,j,k,2) *= rhoinv;);
                });
            } // mfi
        } // lev
    } // not m_advect_momentum
}
