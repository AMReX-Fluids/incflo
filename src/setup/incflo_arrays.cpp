#include <incflo.H>

using namespace amrex;

incflo::LevelData::LevelData (amrex::BoxArray const& ba,
                              amrex::DistributionMapping const& dm,
                              amrex::FabFactory<FArrayBox> const& fact,
                              int ntrac, int ng_state,
                              const std::string& advection_type, bool implicit_diffusion,
                              bool use_tensor_correction, bool advect_tracer)
    : velocity  (ba, dm, AMREX_SPACEDIM, ng_state, MFInfo(), fact),
      velocity_o(ba, dm, AMREX_SPACEDIM, ng_state, MFInfo(), fact),
      velocity_eb(ba, dm, AMREX_SPACEDIM, ng_state, MFInfo(), fact),
      density   (ba, dm, 1             , ng_state, MFInfo(), fact),
      density_eb(ba, dm, 1             , ng_state, MFInfo(), fact),
      density_o (ba, dm, 1             , ng_state, MFInfo(), fact),
      tracer    (ba, dm, ntrac         , ng_state, MFInfo(), fact),
      tracer_eb (ba, dm, ntrac         , ng_state, MFInfo(), fact),
      tracer_o  (ba, dm, ntrac         , ng_state, MFInfo(), fact),
      mac_phi   (ba, dm, 1             , 1       , MFInfo(), fact),
      p_nd      (amrex::convert(ba,IntVect::TheNodeVector()),
                     dm, 1             , 0 , MFInfo(), fact),
      gp        (ba, dm, AMREX_SPACEDIM, 0       , MFInfo(), fact),
      conv_velocity_o(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact),
      conv_density_o (ba, dm, 1             , 0, MFInfo(), fact),
      conv_tracer_o  (ba, dm, ntrac         , 0, MFInfo(), fact),
      energy         (ba, dm, 1, ng_state, MFInfo(), fact),
      energy_o       (ba, dm, 1, ng_state, MFInfo(), fact),
      temp           (ba, dm, 1, ng_state, MFInfo(), fact),
      temp_o         (ba, dm, 1, ng_state, MFInfo(), fact),
      vol_heat_cap   (ba, dm, 1, ng_state, MFInfo(), fact),
      ther_conduct   (ba, dm, 1, ng_state, MFInfo(), fact),
      one            (ba, dm, 1, ng_state, MFInfo(), fact),
      conv_energy_o  (ba, dm, 1             , 0, MFInfo(), fact)
{
    if (advection_type != "MOL") {
        divtau_o.define(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
        // TODO: Add a flag here for cryo-plunging
        // if (advect_tracer) {
            laps_o.define(ba, dm, ntrac, 0, MFInfo(), fact);
        // }
    } else {
        conv_velocity.define(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
        conv_density.define (ba, dm, 1             , 0, MFInfo(), fact);
        conv_tracer.define  (ba, dm, ntrac         , 0, MFInfo(), fact);
        conv_energy.define  (ba, dm, 1             , 0, MFInfo(), fact);

        if (!implicit_diffusion || use_tensor_correction)
        {
            divtau.define  (ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
            divtau_o.define(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
        }
        if (!implicit_diffusion && advect_tracer)
        {
            laps.define  (ba, dm, ntrac, 0, MFInfo(), fact);
            laps_o.define(ba, dm, ntrac, 0, MFInfo(), fact);
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
