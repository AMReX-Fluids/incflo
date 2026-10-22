#include <incflo.H>

using namespace amrex;

//
// Source term S in the divergence constraint div(u) = S.
//
// S^n is the S^(n+1) computed at the end of the previous step: the projection does not
// change the scalars, so there is nothing to recompute at the start of a step.  The MAC
// projection and the Godunov edge-state prediction use S^(n+1/2) = S^n + dt/2 dS/dt,
// where dS/dt is lagged from the previous step, and the nodal projection uses S^(n+1).
// This follows IAMR.
//

void
incflo::compute_divu_source (Vector<MultiFab      *> const& divu_source,
                             Vector<MultiFab const*> const& /*density*/,
                             Vector<MultiFab const*> const& /*tracer*/,
                             Vector<MultiFab const*> const& /*temperature*/,
                             Real time)
{
    BL_PROFILE("incflo::compute_divu_source()");

    for (int lev = 0; lev <= finest_level; ++lev) {
        prob_divu_source(lev, *divu_source[lev], time);
    }
}

// Compute S from the new-time state into divu_source.  The coarse data under the fine
// levels is replaced by the average of the fine data, as for the state.
void
incflo::update_divu_source (Real time)
{
    if (!m_use_divu_source) { return; }

    // Fill one ghost cell of the state, so that compute_divu_source can use a stencil
    for (int lev = 0; lev <= finest_level; ++lev) {
        auto& ld = *m_leveldata[lev];
        fillpatch_density(lev, m_t_new[lev], ld.density, 1);
        if (m_advect_tracer) {
            fillpatch_tracer(lev, m_t_new[lev], ld.tracer, 1);
        }
        if (m_use_temperature) {
            fillpatch_temperature(lev, m_t_new[lev], ld.temperature, 1);
        }
    }

    compute_divu_source(get_divu_source_new(), get_density_new_const(),
                        get_tracer_new_const(), get_temperature_new_const(), time);

    for (int lev = finest_level-1; lev >= 0; --lev) {
        amrex::average_down(m_leveldata[lev+1]->divu_source, m_leveldata[lev]->divu_source,
                            0, 1, refRatio(lev));
    }
}

void
incflo::compute_divu_source_nph ()
{
    if (!m_use_divu_source) { return; }

    Real const half_dt = Real(0.5) * m_dt;
    for (int lev = 0; lev <= finest_level; ++lev) {
        auto& ld = *m_leveldata[lev];
        MultiFab::LinComb(ld.divu_source_nph, Real(1.0), ld.divu_source_o, 0,
                          half_dt, ld.dsdt, 0, 0, 1, 0);
    }

    // The edge-state prediction reads S in ghost cells, including those at the
    // coarse/fine boundary
    Real const time_nph = m_cur_time + half_dt;
    auto divu_source_nph = get_divu_source_nph();
    for (int lev = 0; lev <= finest_level; ++lev) {
        fillpatch_divu_source(lev, time_nph, *divu_source_nph[lev], divu_source_nph,
                              nghost_divu_source_nph());
    }
}

void
incflo::compute_dsdt ()
{
    if (!m_use_divu_source) { return; }

    Real const dtinv = Real(1.0) / m_dt;
    for (int lev = 0; lev <= finest_level; ++lev) {
        auto& ld = *m_leveldata[lev];
        MultiFab::LinComb(ld.dsdt, dtinv, ld.divu_source, 0,
                          -dtinv, ld.divu_source_o, 0, 0, 1, 0);
    }
}

//
// In a closed domain (only walls and periodic boundaries) the projections are singular,
// and MLMG silently removes the mean of the right-hand side, so the velocity satisfies
// div(u) = S - mean(S).  With a uniform 1/(gamma p0) that is the closed-chamber low Mach
// number correction, S - dp0/dt / (gamma p0), but incflo neither evolves p0 nor adds the
// compression heating dp0/dt / (rho cp) to the temperature equation.  Warn once.
//
void
incflo::check_divu_source_solvability ()
{
    if (!m_use_divu_source) { return; }

    bool closed = true;
    for (OrientationIter oit; oit; ++oit) {
        auto const bct = m_bc_type[oit()];
        if (bct != BC::periodic && bct != BC::slip_wall && bct != BC::no_slip_wall) {
            closed = false;
        }
    }
    if (!closed) { return; }
    if (m_divu_solvability_warned && m_verbose <= 0) { return; }

    Vector<MultiFab> abs_S(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        abs_S[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), Factory(lev));
        MultiFab::Copy(abs_S[lev], m_leveldata[lev]->divu_source, 0, 0, 1, 0);
        abs_S[lev].abs(0, 1);
    }
    Real const sum_S     = amrex::volumeWeightedSum(GetVecOfConstPtrs(get_divu_source_new()),
                                                    0, Geom(), refRatio());
    Real const sum_abs_S = amrex::volumeWeightedSum(GetVecOfConstPtrs(abs_S), 0, Geom(),
                                                    refRatio());

    if (m_verbose > 0) {
        amrex::Print() << "Closed domain: integral of S = " << sum_S << "\n";
    }

    constexpr Real tol = Real(1.e-8);
    if (!m_divu_solvability_warned && std::abs(sum_S) > tol * sum_abs_S) {
        m_divu_solvability_warned = true;
        amrex::Print() << "WARNING: incflo.use_divu_source: the domain is closed but the "
                          "integral of S is not zero.  The projections remove the mean of S, "
                          "which is a uniform change of the thermodynamic pressure; incflo "
                          "does not evolve that pressure or add the matching compression "
                          "heating to the temperature equation.\n";
    }
}
