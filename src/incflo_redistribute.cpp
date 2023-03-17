#ifdef AMREX_USE_EB

#include <incflo.H>
#include <hydro_redistribution.H>

using namespace amrex;

void
incflo::redistribute_term ( MultiFab& result,
			    MultiFab& result_tmp, // Saves doing a MF::copy. does this matter???
			    MultiFab const& state,
			    BCRec const* bc, // this is bc for the state (needed for SRD slopes)
			    int lev)
{
    // ************************************************************************
    // Redistribute result_tmp and pass out result
    // ************************************************************************
    AMREX_ASSERT(result.nComp() == state.nComp());

    result_tmp.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
	redistribute_term(mfi, result, result_tmp, state, bc, lev);
    }
}

void
incflo::redistribute_term ( MFIter const& mfi,
			    MultiFab& result,
			    MultiFab& result_tmp,
			    MultiFab const& state,
			    BCRec const* bc, // this is bc for the state (needed for SRD slopes)
			    int lev)
{
    AMREX_ASSERT(result.nComp() == state.nComp());

    Box const& bx = mfi.tilebox();

    EBFArrayBoxFactory const& ebfact = EBFactory(lev);
    EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
    Array4<EBCellFlag const> const& flag = flagfab.const_array();

    bool regular = (flagfab.getType(amrex::grow(bx,4)) == FabType::regular);
    bool covered = (flagfab.getType(bx) == FabType::covered);

    Array4<Real      > out = result.array(mfi);
    Array4<Real      > in  = result_tmp.array(mfi);
    int ncomp = result.nComp();

    if (!regular && !covered)
    {
	auto const& vfrac = ebfact.getVolFrac().const_array(mfi);
	auto const& ccc   = ebfact.getCentroid().const_array(mfi);
	AMREX_D_TERM(auto const& apx = ebfact.getAreaFrac()[0]->const_array(mfi);,
		     auto const& apy = ebfact.getAreaFrac()[1]->const_array(mfi);,
		     auto const& apz = ebfact.getAreaFrac()[2]->const_array(mfi););
	AMREX_D_TERM(auto const& fcx = ebfact.getFaceCent()[0]->const_array(mfi);,
		     auto const& fcy = ebfact.getFaceCent()[1]->const_array(mfi);,
		     auto const& fcz = ebfact.getFaceCent()[2]->const_array(mfi););

	Box gbx = bx;

	if (m_redistribution_type == "StateRedist") {
	    gbx.grow(3);
	} else if (m_redistribution_type == "FluxRedist") {
	    gbx.grow(2);
	}

	FArrayBox scratch_fab(gbx,ncomp);
	Array4<Real> scratch = scratch_fab.array();
	Elixir eli_scratch = scratch_fab.elixir();

	// This is scratch space if calling StateRedistribute
	//  but is used as the weights (here set to 1) if calling
	//  FluxRedistribute
	amrex::ParallelFor(Box(scratch),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
	    scratch(i,j,k) = 1.;
	});
	
	Array4<Real const> state_arr = state.const_array(mfi);
	// State redist acts on a state. What would that be for the diffusive term??
	Redistribution::Apply(bx, ncomp, out, in, state_arr,
			      scratch, flag,
			      AMREX_D_DECL(apx, apy, apz), vfrac,
			      AMREX_D_DECL(fcx, fcy, fcz), ccc,
			      bc, geom[lev], m_dt, m_redistribution_type);
    }
    else
    {
	amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
	    out(i,j,k,n) = in(i,j,k,n);
	});
    }
}
#endif
