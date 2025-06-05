#include <incflo.H>

using namespace amrex;

incflo::LevelData::LevelData (amrex::BoxArray const& ba,
                              amrex::DistributionMapping const& dm,
                              amrex::FabFactory<FArrayBox> const& fact,
                              incflo* my_incflo)
    : velocity    (ba, dm, AMREX_SPACEDIM, my_incflo->nghost_state(), MFInfo(), fact),
      velocity_o  (ba, dm, AMREX_SPACEDIM, my_incflo->nghost_state(), MFInfo(), fact),

      density     (ba, dm, 1             , my_incflo->nghost_state(), MFInfo(), fact),
      density_o   (ba, dm, 1             , my_incflo->nghost_state(), MFInfo(), fact),
      density_nph (ba, dm, 1             , my_incflo->nghost_state(), MFInfo(), fact),

      tracer    (ba, dm, my_incflo->m_ntrac, my_incflo->nghost_state(), MFInfo(), fact),
      tracer_o  (ba, dm, my_incflo->m_ntrac, my_incflo->nghost_state(), MFInfo(), fact),

      mac_phi   (ba, dm, 1             , 1       , MFInfo(), fact),
      gp        (ba, dm, AMREX_SPACEDIM, 0 , MFInfo(), fact),

      conv_velocity_o (ba, dm, AMREX_SPACEDIM    , 0, MFInfo(), fact),
      conv_density_o  (ba, dm, 1                 , 0, MFInfo(), fact),
      conv_tracer_o   (ba, dm, my_incflo->m_ntrac, 0, MFInfo(), fact)
{
    if (my_incflo->m_use_cc_proj) {
        p_cc.define(ba                                  , dm, 1, 1, MFInfo(), fact);
    } else {
        p_nd.define(convert(ba,IntVect::TheNodeVector()), dm, 1, 0, MFInfo(), fact);
    }
    if (my_incflo->m_use_temperature) {
        temperature.define   (ba, dm, 1, my_incflo->nghost_state(), MFInfo(), fact);
        temperature_o.define (ba, dm, 1, my_incflo->nghost_state(), MFInfo(), fact);

        conv_temperature_o.define(ba, dm, 1, 0, MFInfo(), fact);
    }
#ifdef AMREX_USE_EB
    if (my_incflo->hasEBFlow()) {
        velocity_eb.define(ba, dm, AMREX_SPACEDIM, my_incflo->nghost_state(), MFInfo(), fact);
        density_eb.define (ba, dm, 1             , my_incflo->nghost_state(), MFInfo(), fact);
    }
    // Allow for Dirichlet BC on EB even if there's no flow through the EB
    if (my_incflo->m_advect_tracer && !my_incflo->m_eb_flow.tracer.empty()) {
        tracer_eb.define  (ba, dm, my_incflo->m_ntrac, my_incflo->nghost_state(), MFInfo(), fact);
    }
    if (my_incflo->m_use_temperature && !my_incflo->m_eb_flow.temperature.empty()) {
        temperature_eb.define(ba, dm, 1, my_incflo->nghost_state(), MFInfo(), fact);
    }
#endif
    if (my_incflo->m_advection_type != "MOL") {
        divtau_o.define(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
        if (my_incflo->m_advect_tracer) {
            laps_o.define(ba, dm, my_incflo->m_ntrac, 0, MFInfo(), fact);
        }
        if (my_incflo->m_use_temperature) {
            laps_tem_o.define(ba, dm, 1, 0, MFInfo(), fact);
        }
    } else {
        conv_velocity.define(ba, dm, AMREX_SPACEDIM   , 0, MFInfo(), fact);
        conv_density.define (ba, dm, 1                , 0, MFInfo(), fact);
        conv_tracer.define (ba, dm, my_incflo->m_ntrac, 0, MFInfo(), fact);

        if (my_incflo->m_use_temperature) {
            conv_temperature.define(ba, dm, 1, 0, MFInfo(), fact);
        }

        bool implicit_diffusion = my_incflo->m_diff_type == DiffusionType::Implicit;
        if (!implicit_diffusion || my_incflo->use_tensor_correction)
        {
            divtau.define  (ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
            divtau_o.define(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
        }
        if (!implicit_diffusion)
        {
            if ( my_incflo->m_advect_tracer) {
                laps.define  (ba, dm, my_incflo->m_ntrac, 0, MFInfo(), fact);
                laps_o.define(ba, dm, my_incflo->m_ntrac, 0, MFInfo(), fact);
            }
            if (my_incflo->m_use_temperature) {
                laps_tem.define  (ba, dm, 1, 0, MFInfo(), fact);
                laps_tem_o.define(ba, dm, 1, 0, MFInfo(), fact);
            }
        }
    }
}

// Resize all arrays when instance of incflo class is constructed.
// This is only done at the very start of the simulation.
void incflo::ResizeArrays ()
{
    // Time holders for fillpatch stuff
    m_t_new.resize(max_level + 1);
    m_t_old.resize(max_level + 1);

    m_leveldata.resize(max_level+1);

    m_factory.resize(max_level+1);
}
