#include <incflo.H>

using namespace amrex;

incflo::LevelData::LevelData (amrex::BoxArray const& ba,
                              amrex::DistributionMapping const& dm,
                              amrex::FabFactory<FArrayBox> const& fact,
                              incflo* my_incflo)
    : velocity    (ba, dm, AMREX_SPACEDIM, my_incflo->nghost_state(), MFInfo(), fact),
      velocity_o  (ba, dm, AMREX_SPACEDIM, my_incflo->nghost_state(), MFInfo(), fact),
      velocity_eb (ba, dm, AMREX_SPACEDIM, my_incflo->nghost_state(), MFInfo(), fact),

      density     (ba, dm, 1             , my_incflo->nghost_state(), MFInfo(), fact),
      density_eb  (ba, dm, 1             , my_incflo->nghost_state(), MFInfo(), fact),
      density_o   (ba, dm, 1             , my_incflo->nghost_state(), MFInfo(), fact),
      density_nph (ba, dm, 1             , my_incflo->nghost_state(), MFInfo(), fact),

      tracer    (ba, dm, my_incflo->m_ntrac, my_incflo->nghost_state(), MFInfo(), fact),
      tracer_eb (ba, dm, my_incflo->m_ntrac, my_incflo->nghost_state(), MFInfo(), fact),
      tracer_o  (ba, dm, my_incflo->m_ntrac, my_incflo->nghost_state(), MFInfo(), fact),

      mac_phi   (ba, dm, 1             , 1       , MFInfo(), fact),
      p_cc      (ba, dm, 1             , 1       , MFInfo(), fact),
      p_nd      (amrex::convert(ba,IntVect::TheNodeVector()),
                 dm, 1             , 0 , MFInfo(), fact),
      gp        (ba, dm, AMREX_SPACEDIM, 0 , MFInfo(), fact),

      conv_velocity_o (ba, dm, AMREX_SPACEDIM    , 0, MFInfo(), fact),
      conv_density_o  (ba, dm, 1                 , 0, MFInfo(), fact),
      conv_tracer_o   (ba, dm, my_incflo->m_ntrac, 0, MFInfo(), fact)
{
    if (my_incflo->m_advection_type != "MOL") {
        divtau_o.define(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
        if (my_incflo->m_advect_tracer) {
            laps_o.define(ba, dm, my_incflo->m_ntrac, 0, MFInfo(), fact);
        }
    } else {
        conv_velocity.define(ba, dm, AMREX_SPACEDIM   , 0, MFInfo(), fact);
        conv_density.define (ba, dm, 1                , 0, MFInfo(), fact);
        conv_tracer.define (ba, dm, my_incflo->m_ntrac, 0, MFInfo(), fact);

        bool implicit_diffusion = my_incflo->m_diff_type == DiffusionType::Implicit;
        if (!implicit_diffusion || my_incflo->use_tensor_correction)
        {
            divtau.define  (ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
            divtau_o.define(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
        }
        if (!implicit_diffusion && my_incflo->m_advect_tracer)
        {
            laps.define  (ba, dm, my_incflo->m_ntrac, 0, MFInfo(), fact);
            laps_o.define(ba, dm, my_incflo->m_ntrac, 0, MFInfo(), fact);
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
