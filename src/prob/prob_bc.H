#ifndef PROB_BC_H_
#define PROB_BC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

struct IncfloVelFill
{
    int probtype;
    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr IncfloVelFill (int a_probtype,
                             amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>,
                             AMREX_SPACEDIM*2> const& a_bcv)
        : probtype(a_probtype), bcv(a_bcv) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& vel,
                     const int dcomp, const int num_comp,
                     amrex::GeometryData const& geom, const amrex::Real /*time*/,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        // do something for external Dirichlet (amrex::BCType::ext_dir)
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

            if (bc.lo(0) == amrex::BCType::ext_dir && i < domain_box.smallEnd(0))
            {
                vel(i,j,k,dcomp+nc) = bcv[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)][orig_comp+nc];

                if (orig_comp+nc == 0) // x velocity
                {
                    if (31 == probtype)
                    {
                        amrex::Real y = amrex::Real(j+0.5)*(amrex::Real(1.0)/domain_box.length(1));
                        vel(i,j,k,dcomp+nc) *= amrex::Real(6.) * y * (amrex::Real(1.0)-y);
                    }
#if (AMREX_SPACEDIM == 3)
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
#endif
                }
            }
            else if (bc.hi(0) == amrex::BCType::ext_dir && i > domain_box.bigEnd(0))
            {
                vel(i,j,k,dcomp+nc) = bcv[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)][orig_comp+nc];
            }

            if (bc.lo(1) == amrex::BCType::ext_dir && j < domain_box.smallEnd(1))
            {
                vel(i,j,k,dcomp+nc) = bcv[amrex::Orientation(amrex::Direction::y,amrex::Orientation::low)][orig_comp+nc];

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
            else if (bc.hi(1) == amrex::BCType::ext_dir && j > domain_box.bigEnd(1))
            {
                if (16 == probtype)
                {
                    if (orig_comp+nc == 0) { // x velocity
                        amrex::Real x = amrex::Real(i+0.5)*(amrex::Real(1.0)/domain_box.length(0));
                        vel(i,j,k,dcomp+nc) = 16.0 * (x*x*x*x - 2.0 * x*x*x + x*x);
                    }
                } else {
                   vel(i,j,k,dcomp+nc) = bcv[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)][orig_comp+nc];
                }
            }

#if (AMREX_SPACEDIM == 3)
            if (bc.lo(2) == amrex::BCType::ext_dir && k < domain_box.smallEnd(2))
            {
                vel(i,j,k,dcomp+nc) = bcv[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)][orig_comp+nc];

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
            else if (bc.hi(2) == amrex::BCType::ext_dir && k > domain_box.bigEnd(2))
            {
                vel(i,j,k,dcomp+nc) = bcv[amrex::Orientation(amrex::Direction::z,amrex::Orientation::high)][orig_comp+nc];
            }
#endif
        }
    }
};

struct IncfloDenFill
{
    int probtype;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr IncfloDenFill (int a_probtype, amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2> const& a_bcv)
        : probtype(a_probtype), bcv(a_bcv) {}

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

        if (bc.lo(0) == amrex::BCType::ext_dir && i < domain_box.smallEnd(0))
        {
            rho(i,j,k) = bcv[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)];
        }
        else if (bc.hi(0) == amrex::BCType::ext_dir && i > domain_box.bigEnd(0))
        {
            rho(i,j,k) = bcv[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)];
        }

        if (bc.lo(1) == amrex::BCType::ext_dir && j < domain_box.smallEnd(1))
        {
            rho(i,j,k) = bcv[amrex::Orientation(amrex::Direction::y,amrex::Orientation::low)];
        }
        else if (bc.hi(1) == amrex::BCType::ext_dir && j > domain_box.bigEnd(1))
        {
            rho(i,j,k) = bcv[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)];
        }
#if (AMREX_SPACEDIM == 3)
        if (bc.lo(2) == amrex::BCType::ext_dir && k < domain_box.smallEnd(2))
        {
            rho(i,j,k) = bcv[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)];
        }
        else if (bc.hi(2) == amrex::BCType::ext_dir && k > domain_box.bigEnd(2))
        {
            rho(i,j,k) = bcv[amrex::Orientation(amrex::Direction::z,amrex::Orientation::high)];
        }
#endif
    }
};

struct IncfloTracFill
{
    int probtype, ntrac;
    amrex::GpuArray<amrex::Real const*, AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr IncfloTracFill (int a_probtype, int a_ntrac,
                              amrex::GpuArray<amrex::Real const*, AMREX_SPACEDIM*2> const& a_bcv)
        : probtype(a_probtype), ntrac(a_ntrac), bcv(a_bcv) {}

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

            if (bc.lo(0) == amrex::BCType::ext_dir && i < domain_box.smallEnd(0))
            {
                tracer(i,j,k,n) = bcv[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)][n];
            }
            else if (bc.hi(0) == amrex::BCType::ext_dir && i > domain_box.bigEnd(0))
            {
                tracer(i,j,k,n) = bcv[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)][n];
            }

            if (bc.lo(1) == amrex::BCType::ext_dir && j < domain_box.smallEnd(1))
            {
                tracer(i,j,k,n) = bcv[amrex::Orientation(amrex::Direction::y,amrex::Orientation::low)][n];
            }
            else if (bc.hi(1) == amrex::BCType::ext_dir && j > domain_box.bigEnd(1))
            {
                tracer(i,j,k,n) = bcv[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)][n];
            }
#if (AMREX_SPACEDIM == 3)
            if (bc.lo(2) == amrex::BCType::ext_dir && k < domain_box.smallEnd(2))
            {
                tracer(i,j,k,n) = bcv[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)][n];
            }
            else if (bc.hi(2) == amrex::BCType::ext_dir && k > domain_box.bigEnd(2))
            {
                tracer(i,j,k,n) = bcv[amrex::Orientation(amrex::Direction::z,amrex::Orientation::high)][n];
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
