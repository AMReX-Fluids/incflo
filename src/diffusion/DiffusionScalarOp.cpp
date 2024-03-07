#include <incflo.H>
#include <DiffusionScalarOp.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB_Redistribution.H>
#include <memory>
#endif

using namespace amrex;

DiffusionScalarOp::DiffusionScalarOp (incflo* a_incflo)
    : m_incflo(a_incflo)
{
    readParameters();

    LPInfo info_solve;
    info_solve.setMaxCoarseningLevel(m_mg_max_coarsening_level);
    LPInfo info_apply;
    info_apply.setMaxCoarseningLevel(0);
#ifdef AMREX_USE_EB
    int finest_level = m_incflo->finestLevel();
    if (!m_incflo->EBFactory(0).isAllRegular())
    {
        Vector<EBFArrayBoxFactory const*> ebfact;
        for (int lev = 0; lev <= finest_level; ++lev) {
            ebfact.push_back(&(m_incflo->EBFactory(lev)));
        }

        m_eb_scal_solve_op = std::make_unique<MLEBABecLap>(m_incflo->Geom(0,finest_level),
                                                 m_incflo->boxArray(0,finest_level),
                                                 m_incflo->DistributionMap(0,finest_level),
                                                 info_solve, ebfact);
        m_eb_scal_solve_op->setMaxOrder(m_mg_maxorder);
        m_eb_scal_solve_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low ),
                                        m_incflo->get_diffuse_scalar_bc(Orientation::high));

        if (!m_incflo->useTensorSolve())
        {
            m_eb_vel_solve_op = std::make_unique<MLEBABecLap>(m_incflo->Geom(0,finest_level),
                                                    m_incflo->boxArray(0,finest_level),
                                                    m_incflo->DistributionMap(0,finest_level),
                                                    info_solve, ebfact);
            m_eb_vel_solve_op->setMaxOrder(m_mg_maxorder);

            // We don't call setDomainBC here because we will need to call it separately for each component
        }

        if (m_incflo->need_divtau())
        {
            m_eb_scal_apply_op = std::make_unique<MLEBABecLap>(m_incflo->Geom(0,finest_level),
                                                     m_incflo->boxArray(0,finest_level),
                                                     m_incflo->DistributionMap(0,finest_level),
                                                     info_apply, ebfact);
            m_eb_scal_apply_op->setMaxOrder(m_mg_maxorder);

            m_eb_scal_apply_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
                                            m_incflo->get_diffuse_scalar_bc(Orientation::high));
        }

        if ( (m_incflo->need_divtau() && !m_incflo->useTensorSolve()) ||
              m_incflo->useTensorCorrection() )
        {
            m_eb_vel_apply_op = std::make_unique<MLEBABecLap>(m_incflo->Geom(0,finest_level),
                                                    m_incflo->boxArray(0,finest_level),
                                                    m_incflo->DistributionMap(0,finest_level),
                                                    info_apply, ebfact);
            m_eb_vel_apply_op->setMaxOrder(m_mg_maxorder);

            // We don't call setDomainBC here because we will need to call it separately for each component
        }
    }
    else
#endif
    {
        m_reg_scal_solve_op = std::make_unique<MLABecLaplacian>(m_incflo->Geom(0,m_incflo->finestLevel()),
                                                      m_incflo->boxArray(0,m_incflo->finestLevel()),
                                                      m_incflo->DistributionMap(0,m_incflo->finestLevel()),
                                                      info_solve);
        m_reg_scal_solve_op->setMaxOrder(m_mg_maxorder);
        m_reg_scal_solve_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
                                         m_incflo->get_diffuse_scalar_bc(Orientation::high));

        if (!m_incflo->useTensorSolve())
        {
            m_reg_vel_solve_op = std::make_unique<MLABecLaplacian>(m_incflo->Geom(0,m_incflo->finestLevel()),
                                                         m_incflo->boxArray(0,m_incflo->finestLevel()),
                                                         m_incflo->DistributionMap(0,m_incflo->finestLevel()),
                                                         info_solve);
            m_reg_vel_solve_op->setMaxOrder(m_mg_maxorder);

            // We don't call setDomainBC here because we will need to call it separately for each component
        }
        if (m_incflo->need_divtau()) {
            m_reg_scal_apply_op = std::make_unique<MLABecLaplacian>(m_incflo->Geom(0,m_incflo->finestLevel()),
                                                          m_incflo->boxArray(0,m_incflo->finestLevel()),
                                                          m_incflo->DistributionMap(0,m_incflo->finestLevel()),
                                                          info_apply);
            m_reg_scal_apply_op->setMaxOrder(m_mg_maxorder);
            m_reg_scal_apply_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
                                             m_incflo->get_diffuse_scalar_bc(Orientation::high));
        }

        if ( (m_incflo->need_divtau() && !m_incflo->useTensorSolve()) ||
              m_incflo->useTensorCorrection() )
        {
            m_reg_vel_apply_op = std::make_unique<MLABecLaplacian>(m_incflo->Geom(0,m_incflo->finestLevel()),
                                                         m_incflo->boxArray(0,m_incflo->finestLevel()),
                                                         m_incflo->DistributionMap(0,m_incflo->finestLevel()),
                                                         info_apply);
            m_reg_vel_apply_op->setMaxOrder(m_mg_maxorder);

            // We don't call setDomainBC here because we will need to call it separately for each component
        }
    }
}

void
DiffusionScalarOp::readParameters ()
{
    ParmParse pp("scalar_diffusion");

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
    pp.query("bottom_solver", m_mg_bottom_solver);
    pp.query("bottom_verbose", m_mg_bottom_verbose);

    pp.query("num_pre_smooth", m_num_pre_smooth);
    pp.query("num_post_smooth", m_num_post_smooth);
}

void
DiffusionScalarOp::diffuse_scalar (Vector<MultiFab*> const& tracer,
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
        amrex::Print() << "Diffusing scalars one at a time ..." << std::endl;
    }

    const int finest_level = m_incflo->finestLevel();

    Vector<MultiFab> rhs(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rhs[lev].define(tracer[lev]->boxArray(), tracer[lev]->DistributionMap(), 1, 0);
    }

#ifdef AMREX_USE_EB
    if (m_eb_scal_solve_op)
    {
        m_eb_scal_solve_op->setScalars(1.0, dt);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_eb_scal_solve_op->setACoeffs(lev, *density[lev]);
        }
    }
    else
#endif
    {
        m_reg_scal_solve_op->setScalars(1.0, dt);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_reg_scal_solve_op->setACoeffs(lev, *density[lev]);
        }
    }

    for (int comp = 0; comp < tracer[0]->nComp(); ++comp)
    {
#ifdef AMREX_USE_EB
        if (m_eb_scal_solve_op)
        {
            if ( m_incflo->m_has_mixedBC && comp>0 ) {
                // Must reset scalars (and Acoef, done below) to reuse solver with Robin BC
                m_eb_scal_solve_op->setScalars(1.0, dt);
            }

            for (int lev = 0; lev <= finest_level; ++lev) {
                Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_scalar_eta_to_faces(lev, comp, *eta[lev]);
                m_eb_scal_solve_op->setBCoeffs(lev, GetArrOfConstPtrs(b), MLMG::Location::FaceCentroid);
            }
        }
        else
#endif
        {
            for (int lev = 0; lev <= finest_level; ++lev) {
                Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_scalar_eta_to_faces(lev, comp, *eta[lev]);
                m_reg_scal_solve_op->setBCoeffs(lev, GetArrOfConstPtrs(b));
            }
        }

        Vector<MultiFab> phi;
        for (int lev = 0; lev <= finest_level; ++lev) {
            phi.emplace_back(*tracer[lev], amrex::make_alias, comp, 1);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& rhs_a = rhs[lev].array(mfi);
                Array4<Real const> const& tra_a = tracer[lev]->const_array(mfi,comp);
                Array4<Real const> const& rho_a = density[lev]->const_array(mfi);
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rhs_a(i,j,k) = rho_a(i,j,k) * tra_a(i,j,k);
                });
            }

#ifdef AMREX_USE_EB
            if (m_eb_scal_solve_op) {
                if ( m_incflo->m_has_mixedBC ) {
                    if ( comp>0 ) {
                        // Must reset Acoef to reuse solver with Robin BC
                        m_eb_scal_solve_op->setACoeffs(lev, *density[lev]);
                    }

                    auto const robin = m_incflo->make_robinBC_MFs(lev, &phi[lev]);

                    m_eb_scal_solve_op->setLevelBC(lev, &phi[lev],
                                                   &robin[0], &robin[1], &robin[2]);
                }
                else {
                    m_eb_scal_solve_op->setLevelBC(lev, &phi[lev]);
                }

                // For when we use the stencil for centroid values
                // m_eb_scal_solve_op->setPhiOnCentroid();
            } else
#endif
            {
                m_reg_scal_solve_op->setLevelBC(lev, &phi[lev]);
            }
        }

#ifdef AMREX_USE_EB
        MLMG mlmg(m_eb_scal_solve_op ? static_cast<MLLinOp&>(*m_eb_scal_solve_op) : static_cast<MLLinOp&>(*m_reg_scal_solve_op));
#else
        MLMG mlmg(*m_reg_scal_solve_op);
#endif

        // The default bottom solver is BiCG
        if (m_mg_bottom_solver == "smoother")
        {
            mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
        }
        else if (m_mg_bottom_solver == "hypre")
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

        mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), m_mg_rtol, m_mg_atol);
    }
}

void
DiffusionScalarOp::diffuse_vel_components (Vector<MultiFab*> const& vel,
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
        amrex::Print() << "Diffusing velocity components one at a time ..." << std::endl;
    }

    AMREX_ASSERT(vel[0]->nComp() == AMREX_SPACEDIM);

    const int finest_level = m_incflo->finestLevel();

    Vector<MultiFab> rhs(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rhs[lev].define(vel[lev]->boxArray(), vel[lev]->DistributionMap(), 1, 0);
    }

    for (int comp = 0; comp < vel[0]->nComp(); ++comp)
    {
        int eta_comp = 0;

#ifdef AMREX_USE_EB
        if (m_eb_vel_solve_op)
        {
            // Because the different components may have different boundary conditions, we need to
            // reset these for each solve
            m_eb_vel_solve_op->setDomainBC(m_incflo->get_diffuse_velocity_bc(Orientation::low ,comp),
                                           m_incflo->get_diffuse_velocity_bc(Orientation::high,comp));

            m_eb_vel_solve_op->setScalars(1.0, dt);
            for (int lev = 0; lev <= finest_level; ++lev) {
                m_eb_vel_solve_op->setACoeffs(lev, *density[lev]);

                if (m_incflo->hasEBFlow()) {
                  MultiFab phi(*m_incflo->get_velocity_eb()[lev], amrex::make_alias, comp, 1);
                  m_eb_vel_solve_op->setEBDirichlet(lev, phi, *eta[lev]);
                } else {
                  m_eb_vel_solve_op->setEBHomogDirichlet(lev, *eta[lev]);
                }
            }

            for (int lev = 0; lev <= finest_level; ++lev) {
                Array<MultiFab,AMREX_SPACEDIM>
                    b = m_incflo->average_scalar_eta_to_faces(lev, eta_comp, *eta[lev]);

                m_eb_vel_solve_op->setBCoeffs(lev, GetArrOfConstPtrs(b), MLMG::Location::FaceCentroid);
            }
        }
        else
#endif
        {
            // Because the different components may have different boundary conditions, we need to
            // reset these for each solve
            m_reg_vel_solve_op->setDomainBC(m_incflo->get_diffuse_velocity_bc(Orientation::low ,comp),
                                            m_incflo->get_diffuse_velocity_bc(Orientation::high,comp));

            m_reg_vel_solve_op->setScalars(1.0, dt);
            for (int lev = 0; lev <= finest_level; ++lev) {
                m_reg_vel_solve_op->setACoeffs(lev, *density[lev]);
            }

            for (int lev = 0; lev <= finest_level; ++lev) {
                Array<MultiFab,AMREX_SPACEDIM>
                    b = m_incflo->average_scalar_eta_to_faces(lev, eta_comp, *eta[lev]);
                m_reg_vel_solve_op->setBCoeffs(lev, GetArrOfConstPtrs(b));
            }
        }

        Vector<MultiFab> phi;
        for (int lev = 0; lev <= finest_level; ++lev) {
            vel[lev]->FillBoundary(m_incflo->Geom(lev).periodicity());
            phi.emplace_back(*vel[lev], amrex::make_alias, comp, 1);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& rhs_a = rhs[lev].array(mfi);
                Array4<Real const> const& vel_a = vel[lev]->const_array(mfi,comp);
                Array4<Real const> const& rho_a = density[lev]->const_array(mfi);
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rhs_a(i,j,k) = rho_a(i,j,k) * vel_a(i,j,k);
                });
            }

#ifdef AMREX_USE_EB
            if (m_eb_vel_solve_op) {

                // For when we use the stencil for centroid values
                // m_eb_vel_solve_op->setPhiOnCentroid();

                if ( m_incflo->m_has_mixedBC ) {
                    auto const robin = m_incflo->make_robinBC_MFs(lev, &phi[lev]);

                    m_eb_vel_solve_op->setLevelBC(lev, &phi[lev],
                                                  &robin[0], &robin[1], &robin[2]);
                }
                else {
                    m_eb_vel_solve_op->setLevelBC(lev, &phi[lev]);
                }
            } else
#endif
            {
                m_reg_vel_solve_op->setLevelBC(lev, &phi[lev]);
            }
        }

#ifdef AMREX_USE_EB
        MLMG mlmg(m_eb_vel_solve_op ? static_cast<MLLinOp&>(*m_eb_vel_solve_op) : static_cast<MLLinOp&>(*m_reg_vel_solve_op));
#else
        MLMG mlmg(*m_reg_vel_solve_op);
#endif

        // The default bottom solver is BiCG
        if (m_mg_bottom_solver == "smoother")
        {
            mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
        }
        else if (m_mg_bottom_solver == "hypre")
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

        mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), m_mg_rtol, m_mg_atol);
    }
}

void DiffusionScalarOp::compute_laps (Vector<MultiFab*> const& a_laps,
                                      Vector<MultiFab const*> const& a_scalar,
                                      Vector<MultiFab const*> const& a_density,
                                      Vector<MultiFab const*> const& a_eta)
{
    BL_PROFILE("DiffusionScalarOp::compute_laps");

    int finest_level = m_incflo->finestLevel();

    // FIXME - Why do we need this copy?
    Vector<MultiFab> scalar(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_ASSERT(a_scalar[lev]->nComp() == a_laps[lev]->nComp());
        scalar[lev].define(a_scalar[lev]->boxArray(),
                           a_scalar[lev]->DistributionMap(),
                           m_incflo->m_ntrac, 1, MFInfo(),
                           a_scalar[lev]->Factory());
        MultiFab::Copy(scalar[lev], *a_scalar[lev], 0, 0, m_incflo->m_ntrac, 1);
    }

#ifdef AMREX_USE_EB
    if (m_eb_scal_apply_op)
    {
        Vector<MultiFab> laps_tmp(finest_level+1);
        int tmp_comp = (m_incflo->m_redistribution_type == "StateRedist") ? 3 : 2;
        for (int lev = 0; lev <= finest_level; ++lev) {
            laps_tmp[lev].define(a_laps[lev]->boxArray(),
                                 a_laps[lev]->DistributionMap(),
                                 m_incflo->m_ntrac, tmp_comp, MFInfo(),
                                 a_laps[lev]->Factory());
            laps_tmp[lev].setVal(0.0);
        }


        // We want to return div (mu grad)) phi
        m_eb_scal_apply_op->setScalars(0.0, -1.0);

        // For when we use the stencil for centroid values
        // m_eb_scal_apply_op->setPhiOnCentroid();

        // FIXME? Can we do the solve together now?
        for (int comp = 0; comp < m_incflo->m_ntrac; ++comp) {
            int eta_comp = comp;

            if ( m_incflo->m_has_mixedBC && comp>0 ){
                // Must reset scalars to use solver with Robin BC
                m_eb_scal_apply_op->setScalars(0.0, -1.0);
            }

            Vector<MultiFab> laps_comp;
            Vector<MultiFab> scalar_comp;
            for (int lev = 0; lev <= finest_level; ++lev) {
                laps_comp.emplace_back(laps_tmp[lev],amrex::make_alias,comp,1);
                scalar_comp.emplace_back(scalar[lev],amrex::make_alias,comp,1);

                Array<MultiFab,AMREX_SPACEDIM>
                    b = m_incflo->average_scalar_eta_to_faces(lev, eta_comp, *a_eta[lev]);

                m_eb_scal_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(b), MLMG::Location::FaceCentroid);

                if ( m_incflo->m_has_mixedBC ) {

                    auto const robin = m_incflo->make_robinBC_MFs(lev, &scalar_comp[lev]);

                    m_eb_scal_apply_op->setLevelBC(lev, &scalar_comp[lev],
                                                   &robin[0], &robin[1], &robin[2]);
                }
                else {
                    m_eb_scal_apply_op->setLevelBC(lev, &scalar_comp[lev]);
                }
            }

            MLMG mlmg(*m_eb_scal_apply_op);
            mlmg.apply(GetVecOfPtrs(laps_comp), GetVecOfPtrs(scalar_comp));
        }

        for(int lev = 0; lev <= finest_level; lev++)
        {
            amrex::single_level_redistribute(laps_tmp[lev],
                                             *a_laps[lev], 0, m_incflo->m_ntrac,
                                             m_incflo->Geom(lev));
            // auto const& bc = m_incflo->get_tracer_bcrec_device_ptr();
            // m_incflo->redistribute_term(*a_laps[lev], laps_tmp[lev], *a_scalar[lev],
            //                 bc, lev);
        }
    }
    else
#endif
    {
        // We want to return div (mu grad)) phi
        m_reg_scal_apply_op->setScalars(0.0, -1.0);

        // This should have no effect since the first scalar is 0
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_reg_scal_apply_op->setACoeffs(lev, *a_density[lev]);
        }

        for (int comp = 0; comp < m_incflo->m_ntrac; ++comp) {

            int eta_comp = comp;

            Vector<MultiFab> laps_comp;
            Vector<MultiFab> scalar_comp;
            for (int lev = 0; lev <= finest_level; ++lev) {
                laps_comp.emplace_back(*a_laps[lev],amrex::make_alias,comp,1);
                scalar_comp.emplace_back(scalar[lev],amrex::make_alias,comp,1);
                Array<MultiFab,AMREX_SPACEDIM>
                    b = m_incflo->average_scalar_eta_to_faces(lev, eta_comp, *a_eta[lev]);

                m_reg_scal_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(b));
                m_reg_scal_apply_op->setLevelBC(lev, &scalar_comp[lev]);
            }

            MLMG mlmg(*m_reg_scal_apply_op);
            mlmg.apply(GetVecOfPtrs(laps_comp), GetVecOfPtrs(scalar_comp));
        }
    }
}

void DiffusionScalarOp::compute_divtau (Vector<MultiFab*> const& a_divtau,
                                        Vector<MultiFab const*> const& a_vel,
                                        Vector<MultiFab const*> const& a_density,
                                        Vector<MultiFab const*> const& a_eta)
{
    BL_PROFILE("DiffusionScalarOp::compute_divtau");

    int finest_level = m_incflo->finestLevel();

    AMREX_ASSERT(a_vel[0]->nComp()    == AMREX_SPACEDIM);
    AMREX_ASSERT(a_divtau[0]->nComp() == AMREX_SPACEDIM);

    Vector<MultiFab> vel(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel[lev].define(a_vel[lev]->boxArray(),
                        a_vel[lev]->DistributionMap(),
                        a_vel[lev]->nComp(), 1, MFInfo(),
                        a_vel[lev]->Factory());
        MultiFab::Copy(vel[lev], *a_vel[lev], 0, 0, a_vel[lev]->nComp(), 1);
    }

#ifdef AMREX_USE_EB
    if (m_eb_vel_apply_op)
    {
        Vector<MultiFab> divtau_tmp(finest_level+1);
        int tmp_comp = (m_incflo->m_redistribution_type == "StateRedist") ? 3 : 2;
        for (int lev = 0; lev <= finest_level; ++lev) {
            divtau_tmp[lev].define(a_divtau[lev]->boxArray(),
                                   a_divtau[lev]->DistributionMap(),
                                   a_divtau[lev]->nComp(), tmp_comp, MFInfo(),
                                   a_divtau[lev]->Factory());
            divtau_tmp[lev].setVal(0.0);
        }

        for (int lev = 0; lev <= finest_level; ++lev) {
            m_eb_vel_apply_op->setEBHomogDirichlet(lev, *a_eta[lev]);
        }
        // We want to return div (mu grad)) phi
        m_eb_vel_apply_op->setScalars(0.0, -1.0);

        // For when we use the stencil for centroid values
        // m_eb_vel_apply_op->setPhiOnCentroid();

        int eta_comp = 0;

        for (int comp = 0; comp < a_divtau[0]->nComp(); ++comp)
        {
            Vector<MultiFab> divtau_single;
            Vector<MultiFab>    vel_single;

            if ( m_incflo->m_has_mixedBC && comp>0 ){
                // Must reset scalars to use solver with Robin BC
                m_eb_vel_apply_op->setScalars(0.0, -1.0);
            }

            // Because the different components may have different boundary conditions, we need to
            // reset these for each solve
            m_eb_vel_apply_op->setDomainBC(m_incflo->get_diffuse_velocity_bc(Orientation::low ,comp),
                                           m_incflo->get_diffuse_velocity_bc(Orientation::high,comp));

            for (int lev = 0; lev <= finest_level; ++lev) {
                divtau_single.emplace_back(divtau_tmp[lev],amrex::make_alias,comp,1);
                vel_single.emplace_back(       vel[lev],amrex::make_alias,comp,1);

                if ( m_incflo->m_has_mixedBC ) {
                    auto const robin = m_incflo->make_robinBC_MFs(lev, &vel_single[lev]);

                    m_eb_vel_apply_op->setLevelBC(lev, &vel_single[lev],
                                                  &robin[0], &robin[1], &robin[2]);
                }
                else {
                    m_eb_vel_apply_op->setLevelBC(lev, &vel_single[lev]);
                }

                Array<MultiFab,AMREX_SPACEDIM> b =
                    m_incflo->average_scalar_eta_to_faces(lev, eta_comp, *a_eta[lev]);
                m_eb_vel_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(b), MLMG::Location::FaceCentroid);
            }

            MLMG mlmg(*m_eb_vel_apply_op);

            mlmg.apply(GetVecOfPtrs(divtau_single), GetVecOfPtrs(vel_single));
        }

        for(int lev = 0; lev <= finest_level; lev++)
        {
            amrex::single_level_redistribute(divtau_tmp[lev],
                                             *a_divtau[lev], 0, a_divtau[lev]->nComp(),
                                             m_incflo->Geom(lev));
            // auto const& bc = m_incflo->get_velocity_bcrec_device_ptr();
            // m_incflo->redistribute_term(*a_divtau[lev], divtau_tmp[lev], *a_vel[lev],
            //                 bc, lev);
        }
    }
    else
#endif
    {
        // We want to return div (mu grad)) phi
        m_reg_vel_apply_op->setScalars(0.0, -1.0);

        // This should have no effect since the first scalar is 0
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_reg_vel_apply_op->setACoeffs(lev, *a_density[lev]);
        }

        int eta_comp = 0;
        Vector<MultiFab> divtau_single;
        Vector<MultiFab>    vel_single;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            Array<MultiFab,AMREX_SPACEDIM>
                b = m_incflo->average_scalar_eta_to_faces(lev, eta_comp, *a_eta[lev]);
            m_reg_vel_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(b));
        }

        for (int comp = 0; comp < a_divtau[0]->nComp(); ++comp)
        {
            // Because the different components may have different boundary conditions, we need to
            // reset these for each solve
            m_reg_vel_apply_op->setDomainBC(m_incflo->get_diffuse_velocity_bc(Orientation::low ,comp),
                                        m_incflo->get_diffuse_velocity_bc(Orientation::high,comp));

            for (int lev = 0; lev <= finest_level; ++lev) {
                divtau_single.emplace_back(*a_divtau[lev],amrex::make_alias,comp,1);
                   vel_single.emplace_back(      vel[lev],amrex::make_alias,comp,1);

                m_reg_vel_apply_op->setLevelBC(lev, &vel_single[lev]);
            }

            MLMG mlmg(*m_reg_vel_apply_op);
            mlmg.apply(GetVecOfPtrs(divtau_single), GetVecOfPtrs(vel_single));
        }
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
    } // not advect_momentum
}
