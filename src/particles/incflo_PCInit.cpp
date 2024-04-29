#include <incflo_PC.H>

#ifdef INCFLO_USE_PARTICLES

#include <AMReX_ParmParse.H>

using namespace amrex;

/*! Read inputs from file */
void incflo_PC::readInputs ()
{
    BL_PROFILE("incflo_PC::readInputs");

    ParmParse pp(m_name);

    m_initialization_type = incflo_ParticleInitializations::init_box_uniform;
    pp.query("initial_distribution_type", m_initialization_type);

    if (m_initialization_type == incflo_ParticleInitializations::init_box_uniform)
    {
        Vector<Real> particle_box_lo(AMREX_SPACEDIM);
        Vector<Real> particle_box_hi(AMREX_SPACEDIM);

        // Defaults
        for (int i = 0; i < AMREX_SPACEDIM; i++) { particle_box_lo[i] = Geom(0).ProbLo(i); }
        for (int i = 0; i < AMREX_SPACEDIM; i++) { particle_box_hi[i] = Geom(0).ProbHi(i); }

        pp.queryAdd("particle_box_lo", particle_box_lo, AMREX_SPACEDIM);
        AMREX_ASSERT(particle_box_lo.size() == AMREX_SPACEDIM);

        pp.queryAdd("particle_box_hi", particle_box_hi, AMREX_SPACEDIM);
        AMREX_ASSERT(particle_box_hi.size() == AMREX_SPACEDIM);

        m_particle_box.setLo(particle_box_lo);
        m_particle_box.setHi(particle_box_hi);

        // We default to placing the particles randomly within each cell,
        // but can override this for regression testing
        place_randomly_in_cells = true;
        pp.query("place_randomly_in_cells", place_randomly_in_cells);
    }

    m_ppc_init = 1;
    pp.query("initial_particles_per_cell", m_ppc_init);

    m_advect_w_flow = (m_name == incfloParticleNames::tracers ? true : false);
    pp.query("advect_with_flow", m_advect_w_flow);

    return;
}

/*! Initialize particles in domain */
void incflo_PC::InitializeParticles (
#ifdef AMREX_USE_EB
                                      EBFArrayBoxFactory const& ebfact
#endif
                                    )
{
    BL_PROFILE("incflo_PC::initializeParticles");

    if (m_initialization_type == incflo_ParticleInitializations::init_box_uniform)
    {
        initializeParticlesUniformDistributionInBox(m_particle_box
#ifdef AMREX_USE_EB
                                                   ,ebfact
#endif
                                                   );
    } else {
        Print() << "Error: " << m_initialization_type
                << " is not a valid initialization for "
                << m_name << " particle species.\n";
        Error("See error message!");
    }
    return;
}

/*! Uniform distribution: the number of particles per grid cell is specified
 *  by "initial_particles_per_cell", and they are randomly distributed. */
void incflo_PC::initializeParticlesUniformDistributionInBox ( const RealBox& particle_init_domain
#ifdef AMREX_USE_EB
                                                             ,EBFArrayBoxFactory const& ebfact
#endif
                                                            )
{
    BL_PROFILE("incflo_PC::initializeParticlesUniformDistributionInBox");

    const int lev = 0;
    const auto dx = Geom(lev).CellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();

    int particles_per_cell = m_ppc_init;

    iMultiFab num_particles( ParticleBoxArray(lev),
                             ParticleDistributionMap(lev),
                             1, 0 );

    // Default number of particles per cell is 0
    num_particles.setVal(0);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        auto num_particles_arr = num_particles[mfi].array();

#ifdef AMREX_USE_EB
        auto const& vf_arr = ebfact.getVolFrac().const_array(mfi);
#endif

        ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = plo[0] + (i + 0.5)*dx[0];
            Real y = plo[1] + (j + 0.5)*dx[1];
            Real z = plo[2] + (k + 0.5)*dx[2];
#ifdef AMREX_USE_EB
            if (vf_arr(i,j,k) > 0.) {
#endif
            if (i%2 == 0 && j%2 == 1) {
            if (particle_init_domain.contains(RealVect(x,y,z))) {
                num_particles_arr(i,j,k) = particles_per_cell;
            }
            }
#ifdef AMREX_USE_EB
            }
#endif
        });
    }

    iMultiFab offsets( ParticleBoxArray(lev),
                       ParticleDistributionMap(lev),
                       1, 0 );
    offsets.setVal(0);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();

        int np = 0;
        {
            int ncell = num_particles[mfi].numPts();
            const int* in = num_particles[mfi].dataPtr();
            int* out = offsets[mfi].dataPtr();
            np = Scan::PrefixSum<int>( ncell,
                                       [=] AMREX_GPU_DEVICE (int i) -> int { return in[i]; },
                                       [=] AMREX_GPU_DEVICE (int i, int const &x) { out[i] = x; },
                                       Scan::Type::exclusive,
                                       Scan::retSum );
        }
        auto offset_arr = offsets[mfi].array();

        auto& particle_tile = DefineAndReturnParticleTile(lev, mfi);
        particle_tile.resize(np);
        auto aos = &particle_tile.GetArrayOfStructs()[0];
        auto& soa = particle_tile.GetStructOfArrays();
        auto* vx_ptr = soa.GetRealData(incflo_ParticlesRealIdxSoA::vx).data();
        auto* vy_ptr = soa.GetRealData(incflo_ParticlesRealIdxSoA::vy).data();
        auto* vz_ptr = soa.GetRealData(incflo_ParticlesRealIdxSoA::vz).data();
        auto* mass_ptr = soa.GetRealData(incflo_ParticlesRealIdxSoA::mass).data();

        const auto num_particles_arr = num_particles[mfi].array();

        auto my_proc = ParallelDescriptor::MyProc();
        Long pid;
        {
            pid = ParticleType::NextID();
            ParticleType::NextID(pid+np);
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( static_cast<Long>(pid + np) < LastParticleID,
                                          "Error: overflow on particle id numbers!" );

        if (place_randomly_in_cells) {

            ParallelForRNG(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                           const RandomEngine& rnd_engine) noexcept
            {
                int start = offset_arr(i,j,k);
                for (int n = start; n < start+num_particles_arr(i,j,k); n++) {
                    Real r[3] = {Random(rnd_engine), Random(rnd_engine), Random(rnd_engine)};
                    Real v[3] = {0.0, 0.0, 0.0};

                    Real x = plo[0] + (i + r[0])*dx[0];
                    Real y = plo[1] + (j + r[1])*dx[1];
                    Real z = plo[2] + (k + r[2])*dx[2];

                    auto& p = aos[n];
                    p.id()  = pid + n;
                    p.cpu() = my_proc;

                    p.pos(0) = x; p.pos(1) = y; p.pos(2) = z;

                    vx_ptr[n] = v[0]; vy_ptr[n] = v[1]; vz_ptr[n] = v[2];

                    mass_ptr[n] = 1.0e-6;
               }
            });

        } else {

            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int start = offset_arr(i,j,k);
                for (int n = start; n < start+num_particles_arr(i,j,k); n++) {
                    Real r[3] = {0.3, 0.7, 0.25};
                    Real v[3] = {0.0, 0.0, 0.0};

                    Real x = plo[0] + (i + r[0])*dx[0];
                    Real y = plo[1] + (j + r[1])*dx[1];
                    Real z = plo[2] + (k + r[2])*dx[2];

                    auto& p = aos[n];
                    p.id()  = pid + n;
                    p.cpu() = my_proc;

                    p.pos(0) = x; p.pos(1) = y; p.pos(2) = z;

                    vx_ptr[n] = v[0]; vy_ptr[n] = v[1]; vz_ptr[n] = v[2];

                    mass_ptr[n] = 1.0e-6;
               }
            });
        }
    }

    ParmParse pp("cylinder");

    Real cyl_radius;
    pp.get("radius",cyl_radius);

    Array<Real,3> cyl_center;
    pp.get("center",cyl_center);
    Real x_ctr = cyl_center[0];
    Real y_ctr = cyl_center[1];

    // Remove particles that are outside of the cylinder
    for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& ptile = ParticlesAt(lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        auto& soa  = ptile.GetStructOfArrays();
        const int n = aos.numParticles();
        auto *p_pbox = aos().data();

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];

            Real x =  p.pos(0) - x_ctr;
            Real y =  p.pos(1) - y_ctr;

            Real r =  std::sqrt(x*x + y*y);

            if (r > cyl_radius) {
                p.id() = -1;
            }
        });
    }
}

#endif
