#ifdef AMREX_USE_EB

#include <incflo.H>
#include <hydro_redistribution.H>

using namespace amrex;

void
incflo::redistribute_term ( MultiFab& result,
			    MultiFab& temporary, // Saves doing a MF::copy. does this matter???
			    MultiFab const& state,
			    BCRec const* bc, // this is bc for the state (needed for SRD slopes)
			    int lev,
			    MultiFab*& vel_eb)
{
    // ************************************************************************
    // Redistribute result_tmp and pass out result
    // ************************************************************************
    AMREX_ASSERT(result.nComp() == state.nComp());

    temporary.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
	redistribute_term(mfi, result, temporary, state, bc, lev, vel_eb);
    }
}

void
incflo::redistribute_term ( MFIter const& mfi,
			    MultiFab& result,
			    MultiFab& temporary,
			    MultiFab const& state,
			    BCRec const* bc, // this is bc for the state (needed for SRD slopes)
			    int lev,
			    MultiFab*& vel_eb)
{
    Array4<Real      > out       = result.array(mfi);
    Array4<Real      > tmp       = temporary.array(mfi);
    Array4<Real const> state_arr = state.const_array(mfi);
    Array4<Real const> v_eb      = (vel_eb) ? vel_eb->const_array(mfi) : Array4<Real const> {};

    redistribute_term(mfi, out, tmp, state_arr, bc, lev, v_eb);
}

void
incflo::redistribute_term ( MFIter const& mfi,
			    Array4<Real       > const& result,
			    Array4<Real       > const& temporary,
			    Array4<Real const > const& state,
			    BCRec const* bc, // this is bc for the state (needed for SRD slopes)
			    int lev,
			    Array4<Real const > const& vel_eb)
{
    AMREX_ASSERT(result.nComp() == state.nComp());

    Box const& bx = mfi.tilebox();

    EBFArrayBoxFactory const& ebfact = EBFactory(lev);
    EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
    Array4<EBCellFlag const> const& flag = flagfab.const_array();

    bool regular = (flagfab.getType(amrex::grow(bx,4)) == FabType::regular);
    bool covered = (flagfab.getType(bx) == FabType::covered);

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

#ifdef AMREX_USE_MOVING_EB
	EBFArrayBoxFactory const& ebfact_old = OldEBFactory(lev);
	EBCellFlagFab const& flagfab_old         = ebfact_old.getMultiEBCellFlagFab()[mfi];
	Array4<EBCellFlag const> const& flag_old = flagfab_old.const_array();
	auto const& vfrac_old = ebfact_old.getVolFrac().const_array(mfi);
	AMREX_D_TERM(auto const& apx_old = ebfact_old.getAreaFrac()[0]->const_array(mfi);,
		     auto const& apy_old = ebfact_old.getAreaFrac()[1]->const_array(mfi);,
		     auto const& apz_old = ebfact_old.getAreaFrac()[2]->const_array(mfi););
	// For creating the MSRD correction term, so at time n
        Array4<Real const> const& bnorm = ebfact_old.getBndryNormal().const_array(mfi);
        Array4<Real const> const& barea = ebfact_old.getBndryArea().const_array(mfi);


	Redistribution::Apply(bx, ncomp, result, temporary, state,
			      scratch, flag_old, flag,
                              AMREX_D_DECL(apx_old, apy_old, apz_old), vfrac_old,
			      AMREX_D_DECL(apx, apy, apz), vfrac,
			      AMREX_D_DECL(fcx, fcy, fcz), ccc,
			      bc, geom[lev], m_dt, m_redistribution_type,
			      vel_eb, bnorm, barea,
			      Redistribution::defaults::srd_max_order,
			      Redistribution::defaults::target_vol_fraction,
			      Array4<Real const> {});

#else
	// State redist acts on a state. What would that be for the diffusive term??
	Redistribution::Apply(bx, ncomp, result, temporary, state,
			      scratch, flag,
			      AMREX_D_DECL(apx, apy, apz), vfrac,
			      AMREX_D_DECL(fcx, fcy, fcz), ccc,
			      bc, geom[lev], m_dt, m_redistribution_type);
#endif
    }
    else
    {
	amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
	    result(i,j,k,n) = state(i,j,k,n);
	});
    }
}
#endif
