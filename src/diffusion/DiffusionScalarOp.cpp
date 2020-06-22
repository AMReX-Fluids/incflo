#include <incflo.H>
#include <DiffusionScalarOp.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB_utils.H>
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

        m_eb_solve_op.reset(new MLEBABecLap(m_incflo->Geom(0,finest_level),
                                            m_incflo->boxArray(0,finest_level),
                                            m_incflo->DistributionMap(0,finest_level),
                                            info_solve, ebfact));
        m_eb_solve_op->setMaxOrder(m_mg_maxorder);
        // m_eb_solve_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
        //                            m_incflo->get_diffuse_scalar_bc(Orientation::high));

        if (m_incflo->need_divtau()) {
            m_eb_apply_op.reset(new MLEBABecLap(m_incflo->Geom(0,finest_level),
                                                m_incflo->boxArray(0,finest_level),
                                                m_incflo->DistributionMap(0,finest_level),
                                                info_apply, ebfact));
            m_eb_apply_op->setMaxOrder(m_mg_maxorder);
        //     m_eb_apply_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
        //                                m_incflo->get_diffuse_scalar_bc(Orientation::high));
        }
    }
    else
#endif
    {
        m_reg_solve_op.reset(new MLABecLaplacian(m_incflo->Geom(0,m_incflo->finestLevel()),
                                                 m_incflo->boxArray(0,m_incflo->finestLevel()),
                                                 m_incflo->DistributionMap(0,m_incflo->finestLevel()),
                                                 info_solve));
        m_reg_solve_op->setMaxOrder(m_mg_maxorder);
        // m_reg_solve_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
        //                             m_incflo->get_diffuse_scalar_bc(Orientation::high));
        if (m_incflo->need_divtau()) {
            m_reg_apply_op.reset(new MLABecLaplacian(m_incflo->Geom(0,m_incflo->finestLevel()),
                                                     m_incflo->boxArray(0,m_incflo->finestLevel()),
                                                     m_incflo->DistributionMap(0,m_incflo->finestLevel()),
                                                     info_apply));
            m_reg_apply_op->setMaxOrder(m_mg_maxorder);
        //     m_reg_apply_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
        //                                 m_incflo->get_diffuse_scalar_bc(Orientation::high));
        }
    }
}

void
DiffusionScalarOp::readParameters ()
{
    ParmParse pp("scalar_diffusion");

    pp.query("verbose", m_verbose);
    pp.query("mg_verbose", m_mg_verbose);
    pp.query("mg_cg_verbose", m_mg_cg_verbose);
    pp.query("mg_max_iter", m_mg_max_iter);
    pp.query("mg_cg_maxiter", m_mg_cg_maxiter);
    pp.query("mg_max_fmg_iter", m_mg_max_fmg_iter);
    pp.query("mg_max_coarsening_level", m_mg_max_coarsening_level);
    pp.query("mg_maxorder", m_mg_maxorder);
    pp.query("mg_rtol", m_mg_rtol);
    pp.query("mg_atol", m_mg_atol);
    pp.query("bottom_solver_type", m_bottom_solver_type);

    pp.query("num_pre_smooth", m_num_pre_smooth);
    pp.query("num_post_smooth", m_num_post_smooth);
}

void
DiffusionScalarOp::diffuse_scalar (Vector<MultiFab*> const& scalar,
                                   Vector<MultiFab*> const& density,
                                   Vector<MultiFab const*> const& eta,
                                   Real dt, bool is_vel)
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
        if (is_vel)
            amrex::Print() << "Diffusing velocity components one at a time ..." << std::endl;
        else
            amrex::Print() << "Diffusing scalars one at a time ..." << std::endl;
    }

    if (is_vel) AMREX_ASSERT(scalar[lev]->nComp() == AMREX_SPACEDIM);

    const int finest_level = m_incflo->finestLevel();

    Vector<MultiFab> rhs(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rhs[lev].define(scalar[lev]->boxArray(), scalar[lev]->DistributionMap(), 1, 0);
    }

#ifdef AMREX_USE_EB
    if (m_eb_solve_op)
    {
        m_eb_solve_op->setScalars(1.0, dt);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_eb_solve_op->setACoeffs(lev, *density[lev]);
        }
        if (is_vel)
            m_eb_solve_op->setDomainBC(m_incflo->get_diffuse_tensor_bc(Orientation::low),
                                       m_incflo->get_diffuse_tensor_bc(Orientation::high));
        else
            m_eb_solve_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
                                       m_incflo->get_diffuse_scalar_bc(Orientation::high));
    }
    else
#endif
    {
        m_reg_solve_op->setScalars(1.0, dt);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_reg_solve_op->setACoeffs(lev, *density[lev]);
        }
        if (is_vel)
            m_reg_solve_op->setDomainBC(m_incflo->get_diffuse_tensor_bc(Orientation::low),
                                        m_incflo->get_diffuse_tensor_bc(Orientation::high));
        else
            m_reg_solve_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
                                        m_incflo->get_diffuse_scalar_bc(Orientation::high));
    }

    for (int comp = 0; comp < scalar[0]->nComp(); ++comp)
    {
#ifdef AMREX_USE_EB
        if (m_eb_solve_op)
        {
            for (int lev = 0; lev <= finest_level; ++lev) {
                Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_scalar_eta_to_faces(lev, comp, *eta[lev]);
                m_eb_solve_op->setBCoeffs(lev, GetArrOfConstPtrs(b), MLMG::Location::FaceCentroid);
            }
        }
        else
#endif
        {
            for (int lev = 0; lev <= finest_level; ++lev) {
                Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_scalar_eta_to_faces(lev, comp, *eta[lev]);
                m_reg_solve_op->setBCoeffs(lev, GetArrOfConstPtrs(b));
            }
        }

        Vector<MultiFab> phi;
        for (int lev = 0; lev <= finest_level; ++lev) {
            phi.emplace_back(*scalar[lev], amrex::make_alias, comp, 1);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& rhs_a = rhs[lev].array(mfi);
                Array4<Real const> const& tra_a = scalar[lev]->const_array(mfi,comp);
                Array4<Real const> const& rho_a = density[lev]->const_array(mfi);
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rhs_a(i,j,k) = rho_a(i,j,k) * tra_a(i,j,k);
                });
            }

#ifdef AMREX_USE_EB
            if (m_eb_solve_op) {
                m_eb_solve_op->setLevelBC(lev, &phi[lev]);

                // For when we use the stencil for centroid values
                // m_eb_solve_op->setPhiOnCentroid();  
            } else
#endif
            {
                m_reg_solve_op->setLevelBC(lev, &phi[lev]);
            }
        }

#ifdef AMREX_USE_EB
        MLMG mlmg(m_eb_solve_op ? static_cast<MLLinOp&>(*m_eb_solve_op) : static_cast<MLLinOp&>(*m_reg_solve_op));
#else
        MLMG mlmg(*m_reg_solve_op);
#endif

        // The default bottom solver is BiCG
        if (m_bottom_solver_type == "smoother")
        {
            mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
        }
        else if (m_bottom_solver_type == "hypre")
        {
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
        }
        // Maximum iterations for MultiGrid / ConjugateGradients
        mlmg.setMaxIter(m_mg_max_iter);
        mlmg.setMaxFmgIter(m_mg_max_fmg_iter);
        mlmg.setCGMaxIter(m_mg_cg_maxiter);
        
        // Verbosity for MultiGrid / ConjugateGradients
        mlmg.setVerbose(m_mg_verbose);
        mlmg.setCGVerbose(m_mg_cg_verbose);

        mlmg.setPreSmooth(m_num_pre_smooth);
        mlmg.setPostSmooth(m_num_post_smooth);

        mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), m_mg_rtol, m_mg_atol);
    }
}

void DiffusionScalarOp::compute_laps (Vector<MultiFab*> const& a_laps,
                                      Vector<MultiFab const*> const& a_scalar,
                                      Vector<MultiFab const*> const& a_density,
                                      Vector<MultiFab const*> const& a_eta,
                                      bool is_vel)
{
    BL_PROFILE("DiffusionScalarOp::compute_laps");

    int finest_level = m_incflo->finestLevel();

    if (is_vel) AMREX_ASSERT(a_scalar[lev]->nComp() == AMREX_SPACEDIM);

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
    if (m_eb_apply_op)
    {
        Vector<MultiFab> laps_tmp(finest_level+1);
        for (int lev = 0; lev <= finest_level; ++lev) {
            laps_tmp[lev].define(a_laps[lev]->boxArray(),
                                 a_laps[lev]->DistributionMap(),
                                 m_incflo->m_ntrac, 2, MFInfo(),
                                 a_laps[lev]->Factory());
            laps_tmp[lev].setVal(0.0);
        }

        if (is_vel)
            m_eb_apply_op->setDomainBC(m_incflo->get_diffuse_tensor_bc(Orientation::low),
                                       m_incflo->get_diffuse_tensor_bc(Orientation::high));
        else
            m_eb_apply_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
                                       m_incflo->get_diffuse_scalar_bc(Orientation::high));

        // We want to return div (mu grad)) phi
        m_eb_apply_op->setScalars(0.0, -1.0);

        // For when we use the stencil for centroid values
        // m_eb_solve_op->setPhiOnCentroid();  

        // This should have no effect since the first scalar is 0
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_eb_apply_op->setACoeffs(lev, *a_density[lev]);
        }

        for (int comp = 0; comp < m_incflo->m_ntrac; ++comp) {
            Vector<MultiFab> laps_comp;
            Vector<MultiFab> scalar_comp;
            for (int lev = 0; lev <= finest_level; ++lev) {
                laps_comp.emplace_back(laps_tmp[lev],amrex::make_alias,comp,1);
                scalar_comp.emplace_back(scalar[lev],amrex::make_alias,comp,1);
                Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_scalar_eta_to_faces(lev, comp, *a_eta[lev]);
                m_eb_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(b), MLMG::Location::FaceCentroid);
                m_eb_apply_op->setLevelBC(lev, &scalar_comp[lev]);
            }

            MLMG mlmg(*m_eb_apply_op);
            mlmg.apply(GetVecOfPtrs(laps_comp), GetVecOfPtrs(scalar_comp));
        }

        for(int lev = 0; lev <= finest_level; lev++)
        {
            amrex::single_level_redistribute(laps_tmp[lev],
                                             *a_laps[lev], 0, m_incflo->m_ntrac,
                                             m_incflo->Geom(lev));
        }
    }
    else
#endif
    {
        if (is_vel)
            m_reg_apply_op->setDomainBC(m_incflo->get_diffuse_tensor_bc(Orientation::low),
                                        m_incflo->get_diffuse_tensor_bc(Orientation::high));
        else
            m_reg_apply_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
                                        m_incflo->get_diffuse_scalar_bc(Orientation::high));

        // We want to return div (mu grad)) phi
        m_reg_apply_op->setScalars(0.0, -1.0);

        // This should have no effect since the first scalar is 0
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_reg_apply_op->setACoeffs(lev, *a_density[lev]);
        }

        for (int comp = 0; comp < m_incflo->m_ntrac; ++comp) {
            Vector<MultiFab> laps_comp;
            Vector<MultiFab> scalar_comp;
            for (int lev = 0; lev <= finest_level; ++lev) {
                laps_comp.emplace_back(*a_laps[lev],amrex::make_alias,comp,1);
                scalar_comp.emplace_back(scalar[lev],amrex::make_alias,comp,1);
                Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_scalar_eta_to_faces(lev, comp, *a_eta[lev]);
                m_reg_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(b));
                m_reg_apply_op->setLevelBC(lev, &scalar_comp[lev]);
            }

            MLMG mlmg(*m_reg_apply_op);
            mlmg.apply(GetVecOfPtrs(laps_comp), GetVecOfPtrs(scalar_comp));
        }
    }
}
