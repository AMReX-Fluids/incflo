incflo Introduction
=====================

incflo is a new massively parallel code for solving the incompressible Navier-Stokes
equations in the presence of complex geometries.

It is built on top of AMReX, a publicly available software framework designed for building
massively parallel block-structured adaptive mesh refinement (AMR)
applications.

incfo relies on the same fundamental physics as in IAMR but does not support subcycling in time.
A few technical details are below:

-  Fluid velocity, density and tracers are defined at cell centroids; pressure is defined at nodes.

-  There are two possible advection algorithms: a Method-Of-Lines (MOL) approach and a Godunov-method algorithm.
   Both use an intermediate MAC projection for face-centered advection velocities

-  Incompressibility of the fluid is imposed through the use of a projection at the end of the time step

-  The representation of the complex geometry uses the embedded boundary, or cut-cell, approach

-  Hybrid parallelization via MPI+X where X = OpenMP for multicore machines, and CUDA/HIP/DCP++ for CPU/GPU systems

-  Parallel I/O using AMReX native I/O or HDF5.

-  Plotfile format supported by AmrVis, VisIt, ParaView, and yt.

incflo heavily leverages AMReX (see https://amrex-codes.github.io/) which is also supported by 
ECP as part of the AMReX Co-Design Center.

