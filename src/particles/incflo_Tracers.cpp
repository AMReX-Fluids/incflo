#include <string>
#include <incflo.H>
#include <incflo_PC.H>

#ifdef INCFLO_USE_PARTICLES

using namespace amrex;

/*! Read tracer particles parameters */
void incflo::readTracerParticlesParams ()
{
    ParmParse pp("incflo");

    m_use_tracer_particles = 0;

    pp.query(std::string("use_"+incfloParticleNames::tracers).c_str(), m_use_tracer_particles);

    if (m_use_tracer_particles) {
        particleData.addName(incfloParticleNames::tracers);
    }
    return;
}

/*! Initialize tracer particles */
void incflo::initializeTracerParticles ( ParGDBBase* a_gdb
#ifdef AMREX_USE_EB
                                        ,EBFArrayBoxFactory const& ebfact
#endif
                                       )
{
    auto& namelist_unalloc( particleData.getNamesUnalloc() );

    for (auto it = namelist_unalloc.begin(); it != namelist_unalloc.end(); ++it) {

        std::string species_name( *it );

        if (species_name == incfloParticleNames::tracers)
        {
            AMREX_ASSERT(m_use_tracer_particles);
            incflo_PC* pc = new incflo_PC(a_gdb, incfloParticleNames::tracers);
#ifdef AMREX_USE_EB
            pc->InitializeParticles(ebfact);
#else
            pc->InitializeParticles();
#endif
            amrex::Print() << "Initialized " << pc->TotalNumberOfParticles() << " tracer particles.\n";
            particleData.pushBack(incfloParticleNames::tracers, pc);
        }
    }

    if (m_use_tracer_particles) namelist_unalloc.remove( incfloParticleNames::tracers );

    return;
}

/*! Evolve tracers particles for one time step*/
void incflo::evolveTracerParticles (AMREX_D_DECL(Vector<MultiFab const*> const& u_mac,
                                                 Vector<MultiFab const*> const& v_mac,
                                                 Vector<MultiFab const*> const& w_mac))
{
    if (m_use_tracer_particles) {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            particleData[incfloParticleNames::tracers]->EvolveParticles(lev, m_dt,
                                                                        AMREX_D_DECL(u_mac[lev],v_mac[lev],w_mac[lev]));
        }
    }
}
#endif
