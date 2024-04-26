#include <incflo_PC.H>

#ifdef INCFLO_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

/*! Evolve particles for one time step */
void incflo_PC::EvolveParticles ( int                                        a_lev,
                                  Real                                       a_dt_lev,
                                  AMREX_D_DECL(
                                  const MultiFab* a_umac, const MultiFab* a_vmac, const MultiFab* a_wmac))
{
    BL_PROFILE("incflo_PC::EvolveParticles()");

    if (m_advect_w_flow) {
        AdvectWithFlow(a_lev, a_dt_lev, AMREX_D_DECL(a_umac, a_vmac, a_wmac));
    }

    Redistribute();
    return;
}

//
// Uses midpoint method to advance particles using flow velocity.
//
void incflo_PC::AdvectWithFlow (int                                 a_lev,
                                Real                                a_dt,
                                const MultiFab* a_umac, const MultiFab* a_vmac, const MultiFab* a_wmac)
{
    BL_PROFILE("incflo_PC::AdvectWithUmac()");
    AMREX_ASSERT(OK(a_lev, a_lev, a_umac[0].nGrow()-1));
    AMREX_ASSERT(a_lev >= 0 && a_lev < GetParticles().size());

    AMREX_D_TERM(AMREX_ASSERT(a_umac.nGrow() >= 1);,
                 AMREX_ASSERT(a_vmac.nGrow() >= 1);,
                 AMREX_ASSERT(a_wmac.nGrow() >= 1););

    const auto      strttime = amrex::second();
    const Geometry    & geom     = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto dxi = geom.InvCellSizeArray();

    bool is_periodic_in_x = geom.isPeriodic(0);
    bool is_periodic_in_y = geom.isPeriodic(1);
    bool is_periodic_in_z = geom.isPeriodic(2);

    for (int ipass = 0; ipass < 2; ipass++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti)
        {
            int grid    = pti.index();
            auto& ptile = ParticlesAt(a_lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            auto& soa  = ptile.GetStructOfArrays();
            const int n = aos.numParticles();
            auto *p_pbox = aos().data();

            Array<ParticleReal*,AMREX_SPACEDIM> v_ptr;
            v_ptr[0] = soa.GetRealData(incflo_ParticlesRealIdxSoA::vx).data();
            v_ptr[1] = soa.GetRealData(incflo_ParticlesRealIdxSoA::vy).data();
            v_ptr[2] = soa.GetRealData(incflo_ParticlesRealIdxSoA::vz).data();

            const FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*a_umac)[grid]),
                                                                  &((*a_vmac)[grid]),
                                                                  &((*a_wmac)[grid])) };

            //array of these pointers to pass to the GPU
            GpuArray<Array4<const Real>, AMREX_SPACEDIM>
                const umacarr {{AMREX_D_DECL((*fab[0]).array(),
                                             (*fab[1]).array(),
                                             (*fab[2]).array() )}};

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }

                ParticleReal v[AMREX_SPACEDIM];
                mac_interpolate(p, plo, dxi, umacarr, v);

                if (ipass == 0) {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        v_ptr[dim][i] = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5)*a_dt*v[dim]);

                        // Reflect off bottom boundary
                        if (dim == 2 && !is_periodic_in_z && p.pos(dim) < plo[2])
                        {
                            p.pos(dim) = 2.0*plo[dim] - p.pos(dim);
                            v_ptr[dim][i] *= -1.0;
                        }
                        // Reflect off top boundary
                        if (dim == 2 && !is_periodic_in_z && p.pos(dim) > phi[2])
                        {
                            p.pos(dim) = 2.0*phi[dim] - p.pos(dim);
                            v_ptr[dim][i] *= -1.0;
                        }
                    }
                } else {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.pos(dim) = v_ptr[dim][i] + static_cast<ParticleReal>(a_dt*v[dim]);
                        v_ptr[dim][i] = v[dim];

                        // Reflect off bottom boundary
                        if (dim == 2 && !is_periodic_in_z && p.pos(dim) < plo[2])
                        {
                            p.pos(dim) = 2.0*plo[dim] - p.pos(dim);
                            v_ptr[dim][i] *= -1.0;
                        }
                        // Reflect off top boundary
                        if (dim == 2 && !is_periodic_in_z && p.pos(dim) > phi[2])
                        {
                            p.pos(dim)     = 2.0*phi[dim] - p.pos(dim);
                            v_ptr[dim][i] *= -1.0;
                        }
                    }
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
    Real z_ctr = cyl_center[2];

    // Remove particles that are outside of the cylindner
    for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti)
    {
        auto& ptile = ParticlesAt(a_lev, pti);
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

    if (m_verbose > 1)
    {
        auto stoptime = amrex::second() - strttime;

#ifdef AMREX_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
                ParallelReduce::Max(stoptime, ParallelContext::IOProcessorNumberSub(),
                                    ParallelContext::CommunicatorSub());

                Print() << "incflo_PC::AdvectWithFlow() time: " << stoptime << '\n';
#ifdef AMREX_LAZY
        });
#endif
    }
}
#endif
