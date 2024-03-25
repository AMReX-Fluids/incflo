.. role:: cpp(code)
   :language: c++

Gridding and Load Balancing
===========================

incflo leverages AMReX's flexibility when it comes to how to decompose the
computational domain into individual rectangular grids, and how to distribute
those grids to MPI ranks.  There can be grids of different sizes,
more than one grid per MPI rank, and different strategies for distributing the grids to MPI ranks.

AMReX's docmumentation on :ref:`amrex:Chap:ManagingGridHierarchy` contains further details.
See :ref:`amrex:sec:grid_creation` for grids are created, i.e. how the :cpp:`BoxArray` on which
:cpp:`MultiFabs` will be built is defined at each level.
See :ref:`amrex:sec:load_balancing` for the strategies AMReX supports for distributing
grids to MPI ranks, i.e. defining the :cpp:`DistributionMapping` with which
:cpp:`MultiFabs` at that level will be built.

When running on multicore machines with OpenMP, we can also control the distribution of
work by setting the size of grid tiles (by defining :cpp:`fabarray.mfiter_tile_size`).
We can also specify the strategy for assigning tiles to OpenMP threads.
See :ref:`amrex:sec:basics:mfiter:tiling` for more about tiling.
