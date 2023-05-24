#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB2.H>
#endif

static amrex::Real volWgtSum (int lev,
                              const amrex::MultiFab* mf,
                              int comp,
#ifdef AMREX_USE_EB
                              const amrex::EBFArrayBoxFactory& ebfactory,
#endif
                              bool local=false)
{
#ifdef AMREX_USE_EB
  const amrex::MultiFab& volfrac =  ebfactory.getVolFrac();
#endif

  amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
  amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (amrex::MFIter mfi(*mf, false); mfi.isValid(); ++mfi) {

    const amrex::Box& bx = mfi.tilebox();

    amrex::Array4<amrex::Real const> const& mf_arr   = mf->const_array(mfi);
#ifdef AMREX_USE_EB
    amrex::Array4<amrex::Real const> const& vfrc_arr = volfrac.const_array(mfi);

    reduce_op.eval(bx, reduce_data, [comp,mf_arr,vfrc_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
    {
      return mf_arr(i,j,k,comp) * vfrc_arr(i,j,k);
    });
#else
    reduce_op.eval(bx, reduce_data, [comp,mf_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
    {
      return mf_arr(i,j,k,comp);
    });
#endif
  }

  ReduceTuple host_tuple = reduce_data.value(reduce_op);

  amrex::Real sum = amrex::get<0>(host_tuple);

  if (!local)
    amrex::ParallelDescriptor::ReduceRealSum(&sum, 1);

  return sum;
}
