MFiX-Exa Introduction
=====================

MFiX-Exa is a new massively parallel code for computing multiphase
flow in which solid particles interact with the gas surrounding them.

It is built on top of AMReX, a publicly available software framework designed for building
massively parallel block-structured adaptive mesh refinement (AMR)
applications.

MFiX-Exa relies on the same fundamental physics as in MFIX but the spatial and
time discretizations differ.  Specifically,

-  Fluid velocity is defined at cell centers; pressure is defined at nodes.

-  The advection algorithm includes an intermediate MAC projection for face-centered advection velocities

-  Incompressibility of the fluid is imposed through the use of a projection at the
   end of the time step rather than through SIMPLE iterations

-  The representation of the complex geometry includes mesh refinement around the 
   fluid-solid interface, which is represented used the embedded boundary, or
   cut-cell approach

-  Parallelization via flat MPI, OpenMP, hybrid MPI/OpenMP, or MPI/MPI.   

-  Porting of MFiX-Exa to GPUs is in progress.

-  Parallel I/O using AMReX native I/O or HDF5.

-  Plotfile format supported by AmrVis, VisIt, ParaView, and yt.

MFiX-Exa is being developed at NETL and LBNL as part of the U.S. Department of Energy's
Exascale Computing Project (ECP).

MFiX-Exa heavily leverages AMReX (see https://amrex-codes.github.io/) which is also supported by 
ECP as part of the AMReX Co-Design Center.

