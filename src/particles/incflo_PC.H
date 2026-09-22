#ifndef INCFLO_PC_H_
#define INCFLO_PC_H_

#ifdef INCFLO_USE_PARTICLES

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <string>
#include <AMReX_Particles.H>

struct incflo_ParticlesIntIdxAoS
{
    enum {
        ncomps = 0
    };
};

struct incflo_ParticlesRealIdxAoS
{
    enum {
        ncomps = 0
    };
};

struct incflo_ParticlesIntIdxSoA
{
    enum {
        ncomps = 0
    };
};

struct incflo_ParticlesRealIdxSoA
{
    enum {
        AMREX_D_DECL(vx = 0,
                     vy,
                     vz),
        mass,
        ncomps
    };
};

namespace incflo_ParticleInitializations
{
    /* list of particle initializations */
    const std::string init_box_uniform  = "box";
}

namespace incfloParticleNames
{
    const std::string tracers = "tracer_particles";
}

class incflo_PC : public amrex::ParticleContainer<  incflo_ParticlesRealIdxAoS::ncomps,  // AoS real attributes
                                                    incflo_ParticlesIntIdxAoS::ncomps,   // AoS integer attributes
                                                    incflo_ParticlesRealIdxSoA::ncomps,  // SoA real attributes
                                                    incflo_ParticlesIntIdxSoA::ncomps,   // SoA integer attributes
                                                    amrex::DefaultAllocator>
{
    public:

        /*! Constructor */
        incflo_PC ( amrex::ParGDBBase* a_gdb,
                    const std::string& a_name = "particles" )
            : amrex::ParticleContainer< incflo_ParticlesRealIdxAoS::ncomps,   // AoS real attributes
                                        incflo_ParticlesIntIdxAoS::ncomps,    // AoS integer attributes
                                        incflo_ParticlesRealIdxSoA::ncomps,   // SoA real attributes
                                        incflo_ParticlesIntIdxSoA::ncomps,    // SoA integer attributes
                                        amrex::DefaultAllocator> (a_gdb)
        {
            BL_PROFILE("incflo_PC::incflo_PC()");
            m_name = a_name;
            // Qualified so that this is a static call: a derived class is not yet
            // constructed here, so virtual dispatch would not reach its override.
            incflo_PC::readInputs();
        }

        /*! Constructor */
        incflo_PC ( const amrex::Geometry&            a_geom,
                    const amrex::DistributionMapping& a_dmap,
                    const amrex::BoxArray&            a_ba,
                    const std::string&                a_name = "particles" )
            : amrex::ParticleContainer< incflo_ParticlesRealIdxAoS::ncomps,   // AoS real attributes
                                        incflo_ParticlesIntIdxAoS::ncomps,    // AoS real attributes
                                        incflo_ParticlesRealIdxSoA::ncomps,   // SoA real attributes
                                        incflo_ParticlesIntIdxSoA::ncomps,    // SoA integer attributes
                                        amrex::DefaultAllocator> ( a_geom, a_dmap, a_ba )
        {
            BL_PROFILE("incflo_PC::incflo_PC()");
            m_name = a_name;
            // Qualified so that this is a static call: a derived class is not yet
            // constructed here, so virtual dispatch would not reach its override.
            incflo_PC::readInputs();
        }

        /*! Initialize particles in domain */
        virtual void InitializeParticles (
#ifdef AMREX_USE_EB
                                          amrex::EBFArrayBoxFactory const& ebfact
#endif
                                         );

        /*! Evolve particles on one level for one time step; the caller is
         *  responsible for calling Redistribute() once all levels are done */
        virtual void EvolveParticles (int,
                                      amrex::Real,
                                      AMREX_D_DECL(const amrex::MultiFab* a_umac,
                                                   const amrex::MultiFab* a_vmac,
                                                   const amrex::MultiFab* a_wmac));

        /*! Get real-type particle attribute names */
        virtual amrex::Vector<std::string> varNames () const
        {
            BL_PROFILE("incflo_PC::varNames()");
            return {AMREX_D_DECL("xvel","yvel","zvel"),"mass"};
        }

        /*! Get real-type particle attribute names */
        virtual amrex::Vector<std::string> meshPlotVarNames () const
        {
            BL_PROFILE("incflo_PC::varNames()");
            return {"mass_density"};
        }

        /*! Uses midpoint method to advance particles using flow velocity. */
        virtual void AdvectWithFlow ( int,
                                      amrex::Real,
                                      AMREX_D_DECL(const amrex::MultiFab* a_umac,
                                                   const amrex::MultiFab* a_vmac,
                                                   const amrex::MultiFab* a_wmac));

        /*! Compute mass density */
        virtual void massDensity ( amrex::MultiFab&, const int&, const int& a_comp = 0) const;

        /*! Compute mesh variable from particles */
        virtual void computeMeshVar(    const std::string&  a_var_name,
                                        amrex::MultiFab&    a_mf,
                                        const int           a_lev) const
        {
            if (a_var_name == "mass_density") {
                massDensity( a_mf, a_lev );
            } else {
                amrex::Abort("this variable not defined in computeMeshVar");
            }
        }

        /*! Specify if particles should advect with flow */
        inline void setAdvectWithFlow (bool a_flag)
        {
            BL_PROFILE("incflo_PC::setAdvectWithFlow()");
            m_advect_w_flow = a_flag;
        }

        // the following functions should ideally be private or protected, but need to be
        // public due to CUDA extended lambda capture rules

        /*! Default particle initialization */
        void initializeParticlesUniformDistributionInBox ( const amrex::RealBox& particle_init_domain
#ifdef AMREX_USE_EB
                                                          ,amrex::EBFArrayBoxFactory const& ebfact
#endif
                                                          );

    protected:

        bool m_advect_w_flow;               /*!< advect with flow velocity */

        amrex::RealBox m_particle_box;      /*!< box within which to place particles */

        std::string m_name;                 /*!< name of this particle species */

        std::string m_initialization_type;  /*!< initial particle distribution type */
        int m_ppc_init;                     /*!< initial number of particles per cell */

        /*! read inputs from file */
        virtual void readInputs ();

    private:

        bool place_randomly_in_cells; /*!< place particles at random positions? */
};

#endif
#endif
