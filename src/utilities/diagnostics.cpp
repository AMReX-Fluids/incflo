#include <incflo.H>

using namespace amrex;

//
// Print maximum values (useful for tracking evolution)
//
void incflo::PrintMaxValues(Real time_in)
{
    return; // xxxxx TODO

#if 0
    ComputeDivU(time_in);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        amrex::Print() << "Level " << lev << std::endl;
        PrintMaxVel(lev);
        PrintMaxGp(lev);
    }
    amrex::Print() << std::endl;
#endif
}

//
// Print the maximum values of the velocity components and velocity divergence
//
void incflo::PrintMaxVel(int lev)
{
#if 0
    // xxxxx TODO
    amrex::Print() << "max(abs(u/v/w/divu))  = "
                   << Norm(vel , lev, 0, 0) << "  "
		   << Norm(vel , lev, 1, 0) << "  "
                   << Norm(vel , lev, 2, 0) << "  "
                   << Norm(divu, lev, 0, 0) << "  " << std::endl;
    if (ntrac > 0)
    {
       for (int i = 0; i < ntrac; i++)
          amrex::Print() << "max tracer" << i << " = " << Norm(tracer,lev,i,0) << std::endl;;
    }
#endif
}

//
// Print the maximum values of the pressure gradient components and pressure
//
void incflo::PrintMaxGp(int lev)
{
#if 0
    // xxxxx TODO
    amrex::Print() << "max(abs(gpx/gpy/gpz/p))  = "
                   << Norm(gp, lev, 0, 0) << "  "
		   << Norm(gp, lev, 1, 0) << "  "
                   << Norm(gp, lev, 2, 0) << "  "
		   << Norm(p , lev, 0, 0) << "  " << std::endl;
#endif
}
