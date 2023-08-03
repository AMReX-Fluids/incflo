#include <incflo.H>

using namespace amrex;

Array<LinOpBCType,AMREX_SPACEDIM>
incflo::get_nodal_projection_bc (Orientation::Side side) const noexcept
{
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = m_bc_type[Orientation(dir,side)];
            switch (bc)
            {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            {
                r[dir] = LinOpBCType::Dirichlet;
                break;
            }
            case BC::mass_inflow:
	    {
                r[dir] = LinOpBCType::inflow;
		break;
            }
            case BC::slip_wall:
            case BC::no_slip_wall:
            {
                r[dir] = LinOpBCType::Neumann;
                break;
            }
            default:
                amrex::Abort("get_nodal_projection_bc: undefined BC type");
            };
        }
    }
    return r;
}

Array<LinOpBCType,AMREX_SPACEDIM>
incflo::get_mac_projection_bc (Orientation::Side side) const noexcept
{
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = m_bc_type[Orientation(dir,side)];
            switch (bc)
            {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            {
                r[dir] = LinOpBCType::Dirichlet;
                break;
            }
            case BC::mass_inflow:
            case BC::slip_wall:
            case BC::no_slip_wall:
            {
                r[dir] = LinOpBCType::Neumann;
                break;
            }
            default:
                amrex::Abort("get_mac_projection_bc: undefined BC type");
            };
        }
    }
    return r;
}
