Introduction
============

incflo is a massively parallel code for solving the incompressible Navier-Stokes
equations in 2-D or 3-D,  with the option for an embedded boundary (cut cell) representation of
of complex geometries.

It is built on top of AMReX, a publicly available software framework designed for building
massively parallel block-structured adaptive mesh refinement (AMR)
applications.

Another AMReX-based code, IAMR, also solves the variable-density incompressible
Navier-Stokes equations in 2-D or 3-D but is based on a subcycling-in-time approach.

Key software and algorithmic features of IAMR include:

-  Fluid velocity, density and tracers are defined at cell centroids; pressure is defined at nodes.

-  Possible advection algorithms: a Method-Of-Lines (MOL) approach and a Godunov-method algorithm.
   Both use an intermediate MAC projection for face-centered advection velocities

-  Incompressibility of the fluid is imposed through the use of a projection at the end of the time step

-  Implicit or explicit discretization of viscous terms with variable viscosity

-  The representation of the complex geometry uses the embedded boundary, or cut-cell, approach

-  Hybrid parallelization via MPI+X where X = OpenMP for multicore machines, and CUDA/HIP/DCP++ for CPU/GPU systems

-  Parallel I/O using AMReX native I/O or HDF5.

-  Plotfile format supported by AmrVis, VisIt, ParaView, and yt.

The incflo source code can be found at https://github.com/AMReX-Codes/incflo/.
incflo heavily leverages AMReX (see https://amrex-codes.github.io/) which is supported by 
ECP as part of the AMReX Co-Design Center.

