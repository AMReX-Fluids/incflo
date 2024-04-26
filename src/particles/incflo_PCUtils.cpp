#include <AMReX_ParticleInterpolators.H>
#include <incflo_PC.H>

#ifdef INCFLO_USE_PARTICLES

using namespace amrex;

void incflo_PC::massDensity ( MultiFab&  a_mf,
                              const int& a_lev,
                              const int& a_comp ) const
{
    BL_PROFILE("incflo_PC::massDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const Real inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    a_mf.setVal(0.0);

    ParticleToMesh( *this, a_mf, a_lev,
        [=] AMREX_GPU_DEVICE (  const incflo_PC::ParticleTileType::ConstParticleTileDataType& ptd,
                                int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh ( p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE ( const incflo_PC::ParticleType&, int)
                {
                    auto mass = ptd.m_rdata[incflo_ParticlesRealIdxSoA::mass][i];
                    return mass*inv_cell_volume;
                });
        });

    return;
}

#endif
