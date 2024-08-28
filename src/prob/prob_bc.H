#ifndef INCFLO_PROB_BC_H_
#define INCFLO_PROB_BC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

struct IncfloVelFill
{
    int probtype;
    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM*2> bcv_vel;

    AMREX_GPU_HOST
    constexpr IncfloVelFill (int a_probtype,
                             amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>,AMREX_SPACEDIM*2> const& a_bcv)
        : probtype(a_probtype), bcv_vel(a_bcv) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& vel,
                     const int dcomp, const int num_comp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        // do something for external Dirichlet (amrex::BCType::ext_dir)
        // and for the Dirichlet part of mixed BCs - which gets FOEXTRAP so the
        // Neumann part is already filled before arriving here during FillPatch.
        const int i = iv[0];
        const int j = iv[1];
#if (AMREX_SPACEDIM == 3)
        const int k = iv[2];
#else
        const int k = 0;
#endif
        const amrex::Box& domain_box = geom.Domain();

        for (int nc = 0; nc < num_comp; ++nc)
        {
            const amrex::BCRec& bc = bcr[bcomp+nc];
#if (AMREX_SPACEDIM == 3)
            else if (1102 == probtype && j > domain_box.bigEnd(1))
            {
                int direction = 2;
                int half_num_cells = domain_box.length(direction) / 2;
                if (k <= half_num_cells) {
                    // Here we take minus the inflow BC specified in inputs file
                    vel(i,j,k,dcomp+nc) = -bcv_vel[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)][orig_comp+nc];
                }
            }
#endif
            // **********************************************************************************************
            // LOW I
            // **********************************************************************************************
            if (i < domain_box.smallEnd(0))
            {
                if (1101 == probtype)
                {
                    int direction = 1;
                    int half_num_cells = domain_box.length(direction) / 2;
                    if (j > half_num_cells) {
                        // Here is the dirichlet portion
                        vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)][orig_comp+nc];
                    }
                }
                else if ( bc.lo(0) == amrex::BCType::ext_dir )
                {
                    vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)][orig_comp+nc];

                    if (orig_comp+nc == 0) // x velocity
                    {
                        if (31 == probtype)
                        {
                            amrex::Real y = amrex::Real(j+0.5)*(amrex::Real(1.0)/domain_box.length(1));
                            vel(i,j,k,dcomp+nc) *= amrex::Real(6.) * y * (amrex::Real(1.0)-y);
                        }
#if (AMREX_SPACEDIM == 2)
                        else if (42 == probtype)
                        {
                            vel(i,j,k,dcomp+nc) = time;
                        }
#elif (AMREX_SPACEDIM == 3)
                        else if (311 == probtype)
                        {
                            amrex::Real z = amrex::Real(k+0.5)*(amrex::Real(1.0)/domain_box.length(2));
                            vel(i,j,k,dcomp+nc) *= amrex::Real(6.) * z * (amrex::Real(1.0)-z);
                        }
                        else if (41 == probtype)
                        {
                            amrex::Real z = amrex::Real(k+0.5)*(amrex::Real(1.0)/domain_box.length(2));
                            vel(i,j,k,dcomp+nc) =  amrex::Real(0.5) * z;
                        }
                        else if (42 == probtype)
                        {
                            vel(i,j,k,dcomp+nc) = time;
                        }
#endif
                    }
                }
                else if ( bc.lo(0) == amrex::BCType::direction_dependent)
                {
                    if (bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)][0] >= 0.)
                    {
                        // do ext_dir if inflowing
                        vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)][orig_comp+nc];
                    } else {
                        // do extrap if outflowing
                        vel(i,j,k,dcomp+nc) = vel(i+1,j,k,dcomp+nc);
                    }
                }
            } // low i

            // **********************************************************************************************
            // HIGH I
            // **********************************************************************************************
            if (i > domain_box.bigEnd(0))
            {
                if (1101 == probtype)
                {
                    int direction = 1;
                    int half_num_cells = domain_box.length(direction) / 2;
                    if (j <= half_num_cells) {
                        // Here we take minus the inflow BC
                        vel(i,j,k,dcomp+nc) = -bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)][orig_comp+nc];
                    }
                }
                else if ( ( (bc.hi(0) == amrex::BCType::ext_dir) ||
                            (bc.hi(0) == amrex::BCType::direction_dependent &&
                             bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)][0] <= 0.) ) )
                {
                    // do ext_dir if inflowing
                    vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)][orig_comp+nc];
                }

                else if ( (bc.hi(0) == amrex::BCType::direction_dependent) &&
                          (bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)][0] > 0.) )
                {
                    // do extrap if outflowing
                    vel(i,j,k,dcomp+nc) = vel(i-1,j,k,dcomp+nc);
                }
            } // high i

            // **********************************************************************************************
            // LOW J
            // **********************************************************************************************
            if ( j < domain_box.smallEnd(1) )
            {
                if ( bc.lo(1) == amrex::BCType::ext_dir )
                {
                    vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::y,amrex::Orientation::low)][orig_comp+nc];

                    if (orig_comp+nc == 1) // y velocity
                    {
#if (AMREX_SPACEDIM == 3)
                        if (32 == probtype)
                        {
                            amrex::Real z = amrex::Real(k+0.5)*(amrex::Real(1.0)/domain_box.length(2));
                            vel(i,j,k,dcomp+nc) *= amrex::Real(6.) * z * (amrex::Real(1.0)-z);
                        }
                        else
#endif
                        if (322 == probtype)
                        {
                            amrex::Real x = amrex::Real(i+0.5)*(amrex::Real(1.0)/domain_box.length(0));
                            vel(i,j,k,dcomp+nc) *= amrex::Real(6.) * x * (amrex::Real(1.0)-x);
                        }
                    }
                }
                else if ( bc.lo(1) == amrex::BCType::direction_dependent )
                {
                    if (bcv_vel[amrex::Orientation(amrex::Direction::y,amrex::Orientation::low)][0] >= 0.)
                    {
                        // do ext_dir if inflowing
                        vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::y,amrex::Orientation::low)][orig_comp+nc];
                    } else {
                        // do extrap if outflowing
                        vel(i,j,k,dcomp+nc) = vel(i,j+1,k,dcomp+nc);
                    }
                }
            } // low j

            // **********************************************************************************************
            // HIGH J
            // **********************************************************************************************
            if ( j > domain_box.bigEnd(1) )
            {
                if ( bc.hi(1) == amrex::BCType::ext_dir )
                {
                    if (16 == probtype)
                    {
                        if (orig_comp+nc == 0) { // x velocity
                            amrex::Real x = amrex::Real(i+0.5)*(amrex::Real(1.0)/domain_box.length(0));
                            vel(i,j,k,dcomp+nc) = 16.0 * (x*x*x*x - 2.0 * x*x*x + x*x);
                        }
                    } else {
                       vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)][orig_comp+nc];
                    }
                }
                else if ( bc.hi(1) == amrex::BCType::direction_dependent )
                {
                    if (bcv_vel[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)][0] <= 0.)
                    {
                        // do ext_dir if inflowing
                        vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)][orig_comp+nc];
                    } else {
                        // do extrap if outflowing
                        vel(i,j,k,dcomp+nc) = vel(i,j-1,k,dcomp+nc);
                    }
                }
            } // high j

#if (AMREX_SPACEDIM == 3)
            // **********************************************************************************************
            // LOW K
            // **********************************************************************************************
            if ( k < domain_box.smallEnd(2) )
            {
                if ( bc.lo(2) == amrex::BCType::ext_dir )
                {
                    vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)][orig_comp+nc];

                    if (orig_comp+nc == 2) // z velocity
                    {
                        if (33 == probtype)
                        {
                            amrex::Real x = amrex::Real(i+0.5)*(amrex::Real(1.0)/domain_box.length(0));
                            vel(i,j,k,dcomp+nc) *= amrex::Real(6.0) * x * (amrex::Real(1.0)-x);
                        }
                        if (333 == probtype)
                        {
                            amrex::Real y = amrex::Real(j+0.5)*(amrex::Real(1.0)/domain_box.length(1));
                            vel(i,j,k,dcomp+nc) *= amrex::Real(6.0) * y * (amrex::Real(1.0)-y);
                        }
                    }
                }
                else if ( bc.lo(2) == amrex::BCType::direction_dependent )
                {
                    if (bcv_vel[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)][0] >= 0.)
                    {
                        // do ext_dir if inflowing
                        vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)][orig_comp+nc];
                    } else {
                        // do extrap if outflowing
                        vel(i,j,k,dcomp+nc) = vel(i,j,k+1,dcomp+nc);
                    }
                }
            } // low k

            // **********************************************************************************************
            // HIGH K
            // **********************************************************************************************
            if ( k > domain_box.bigEnd(2) )
            {
                if ( bc.hi(2) == amrex::BCType::direction_dependent )
                {
                    if (bcv_vel[amrex::Orientation(amrex::Direction::z,amrex::Orientation::high)][0] <= 0.)
                    {
                        // do ext_dir if inflowing
                        vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::z,amrex::Orientation::high)][orig_comp+nc];
                    } else {
                        // do extrap if outflowing
                        vel(i,j,k,dcomp+nc) = vel(i,j,k-1,dcomp+nc);
                    }
                }
                else if ( bc.hi(2) == amrex::BCType::ext_dir )
                {
                    vel(i,j,k,dcomp+nc) = bcv_vel[amrex::Orientation(amrex::Direction::z,amrex::Orientation::high)][orig_comp+nc];
                }
            } // high k
#endif
        } // nc
    }
};

struct IncfloDenFill
{
    int probtype;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM*2> bcv_den;
    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM*2> bcv_vel;

    AMREX_GPU_HOST
    constexpr IncfloDenFill (int a_probtype, amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2> const& a_bcv_den,
                             amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>,AMREX_SPACEDIM*2> const& a_bcv_vel)
        : probtype(a_probtype), bcv_den(a_bcv_den), bcv_vel(a_bcv_vel) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& rho,
                     const int /*dcomp*/, const int /*numcomp*/,
                     amrex::GeometryData const& geom, const amrex::Real /*time*/,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int /*orig_comp*/) const
    {
        using namespace amrex;

        // do something for external Dirichlet (amrex::BCType::ext_dir)
        const int i = iv[0];
        const int j = iv[1];
#if (AMREX_SPACEDIM == 3)
        const int k = iv[2];
#else
        const int k = 0;
#endif

        const Box& domain_box = geom.Domain();
        const BCRec& bc = bcr[bcomp];

        if (1101 == probtype && i < domain_box.smallEnd(0))
        {
            // this probtype uses mixed BC on x faces
            int direction = 1;
            int half_num_cells = domain_box.length(direction) / 2;
            if (j > half_num_cells) {
                // Here is the dirichlet portion; we take the inflow BC specified in inputs file
                rho(i,j,k) = bcv_den[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)];
            }
        }
        else if (1101 == probtype && i > domain_box.bigEnd(0))
        {
            int direction = 1;
            int half_num_cells = domain_box.length(direction) / 2;
            if (j <= half_num_cells) {
                rho(i,j,k) = bcv_den[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)];
            }
        }
#if (AMREX_SPACEDIM == 3)
        else if (1102 == probtype && j > domain_box.bigEnd(1))
        {
            int direction = 2;
            int half_num_cells = domain_box.length(direction) / 2;
            if (k <= half_num_cells) {
                rho(i,j,k) = bcv_den[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)];
            }
        }
#endif
        else if ( (i < domain_box.smallEnd(0)) &&
                  ( (bc.lo(0) == amrex::BCType::ext_dir) ||
                    (bc.lo(0) == amrex::BCType::direction_dependent &&
                     bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)][0] >= 0.) ) )
        {
            rho(i,j,k) = bcv_den[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)];
        }
        else if ( (i > domain_box.bigEnd(0)) &&
                  ( (bc.hi(0) == amrex::BCType::ext_dir) ||
                    (bc.hi(0) == amrex::BCType::direction_dependent &&
                     bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)][0] <= 0.) ) )
        {
            rho(i,j,k) = bcv_den[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)];
        }

        if ( (j < domain_box.smallEnd(1)) &&
              ( (bc.lo(1) == amrex::BCType::ext_dir) ||
                (bc.lo(1) == amrex::BCType::direction_dependent &&
                 bcv_vel[amrex::Orientation(amrex::Direction::y,amrex::Orientation::low)][1] >= 0.) ) )
        {
            rho(i,j,k) = bcv_den[amrex::Orientation(amrex::Direction::y,amrex::Orientation::low)];
        }
        else if ( (j > domain_box.bigEnd(1)) &&
                  ( (bc.hi(1) == amrex::BCType::ext_dir) ||
                    (bc.hi(1) == amrex::BCType::direction_dependent &&
                     bcv_vel[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)][1] <= 0.) ) )
        {
            rho(i,j,k) = bcv_den[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)];
        }
#if (AMREX_SPACEDIM == 3)
         if ( (k < domain_box.smallEnd(2)) &&
              ( (bc.lo(2) == amrex::BCType::ext_dir) ||
                (bc.lo(2) == amrex::BCType::direction_dependent &&
                 bcv_vel[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)][2] >= 0.) ) )
        {
            rho(i,j,k) = bcv_den[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)];
        }
        else if ( (k > domain_box.bigEnd(2)) &&
                  ( (bc.hi(2) == amrex::BCType::ext_dir) ||
                    (bc.hi(2) == amrex::BCType::direction_dependent &&
                     bcv_vel[amrex::Orientation(amrex::Direction::z,amrex::Orientation::high)][2] <= 0.) ) )
        {
            rho(i,j,k) = bcv_den[amrex::Orientation(amrex::Direction::z,amrex::Orientation::high)];
        }
#endif
    }
};

struct IncfloTracFill
{
    int probtype, ntrac;
    amrex::GpuArray<amrex::Real const*, AMREX_SPACEDIM*2> bcv_tra;
    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM*2> bcv_vel;

    AMREX_GPU_HOST
    constexpr IncfloTracFill (int a_probtype, int a_ntrac,
                              amrex::GpuArray<amrex::Real const*, AMREX_SPACEDIM*2> const& a_bcv_tra,
                              amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>,AMREX_SPACEDIM*2> const& a_bcv_vel)
        : probtype(a_probtype), ntrac(a_ntrac), bcv_tra(a_bcv_tra), bcv_vel(a_bcv_vel) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& tracer,
                     const int /*dcomp*/, const int /*numcomp*/,
                     amrex::GeometryData const& geom, const amrex::Real /*time*/,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int /*orig_comp*/) const
    {
        using namespace amrex;

        // do something for external Dirichlet (amrex::BCType::ext_dir)
        const int i = iv[0];
        const int j = iv[1];
#if (AMREX_SPACEDIM == 3)
        const int k = iv[2];
#else
        const int k = 0;
#endif

        const Box& domain_box = geom.Domain();

        for (int n = 0; n < ntrac; ++n)
        {
            const BCRec& bc = bcr[bcomp+n];

            if (1101 == probtype && i < domain_box.smallEnd(0))
            {
                int direction = 1;
                int half_num_cells = domain_box.length(direction) / 2;
                if (j > half_num_cells) {
                    // Here we take the inflow BC specified in inputs file
                    tracer(i,j,k,n) = bcv_tra[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)][n];
                }
            }
            else if (1101 == probtype && i > domain_box.bigEnd(0))
            {
                int direction = 1;
                int half_num_cells = domain_box.length(direction) / 2;
                if (j <= half_num_cells) {
                    tracer(i,j,k,n) = bcv_tra[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)][n];
                }
            }
#if (AMREX_SPACEDIM == 3)
            else if (1102 == probtype && j > domain_box.bigEnd(1))
            {
                int direction = 2;
                int half_num_cells = domain_box.length(direction) / 2;
                if (k <= half_num_cells) {
                    tracer(i,j,k,n) = bcv_tra[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)][n];
                }
            }
#endif
            else if ( (i < domain_box.smallEnd(0)) &&
                      ( (bc.lo(0) == amrex::BCType::ext_dir) ||
                        (bc.lo(0) == amrex::BCType::direction_dependent &&
                         bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)][0] >= 0.) ) )
            {
                tracer(i,j,k,n) = bcv_tra[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)][n];
            }
            else if ( (i > domain_box.bigEnd(0)) &&
                      ( (bc.hi(0) == amrex::BCType::ext_dir) ||
                        (bc.hi(0) == amrex::BCType::direction_dependent &&
                         bcv_vel[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)][0] <= 0.) ) )
            {
                tracer(i,j,k,n) = bcv_tra[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)][n];
            }

            if ( (j < domain_box.smallEnd(1)) &&
                  ( (bc.lo(1) == amrex::BCType::ext_dir) ||
                    (bc.lo(1) == amrex::BCType::direction_dependent &&
                     bcv_vel[amrex::Orientation(amrex::Direction::y,amrex::Orientation::low)][1] >= 0.) ) )
            {
                tracer(i,j,k,n) = bcv_tra[amrex::Orientation(amrex::Direction::y,amrex::Orientation::low)][n];
            }
            else if ( (j > domain_box.bigEnd(1)) &&
                      ( (bc.hi(1) == amrex::BCType::ext_dir) ||
                        (bc.hi(1) == amrex::BCType::direction_dependent &&
                         bcv_vel[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)][1] <= 0.) ) )
            {
                tracer(i,j,k,n) = bcv_tra[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)][n];
            }
#if (AMREX_SPACEDIM == 3)
            if ( (k < domain_box.smallEnd(2)) &&
                 ( (bc.lo(2) == amrex::BCType::ext_dir) ||
                   (bc.lo(2) == amrex::BCType::direction_dependent &&
                    bcv_vel[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)][2] >= 0.) ) )
            {
                tracer(i,j,k,n) = bcv_tra[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)][n];
            }
            else if ( (k > domain_box.bigEnd(2)) &&
                      ( (bc.hi(2) == amrex::BCType::ext_dir) ||
                        (bc.hi(2) == amrex::BCType::direction_dependent &&
                         bcv_vel[amrex::Orientation(amrex::Direction::z,amrex::Orientation::high)][2] <= 0.) ) )
            {
                tracer(i,j,k,n) = bcv_tra[amrex::Orientation(amrex::Direction::z,amrex::Orientation::high)][n];
            }
#endif
        }
    }
};

struct IncfloForFill
{
    int probtype;

    AMREX_GPU_HOST
    constexpr IncfloForFill (int a_probtype) : probtype(a_probtype) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& /*iv*/, amrex::Array4<amrex::Real> const& /*vel*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     amrex::GeometryData const& /*geom*/, const amrex::Real /*time*/,
                     const amrex::BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
    {}
};

#endif
