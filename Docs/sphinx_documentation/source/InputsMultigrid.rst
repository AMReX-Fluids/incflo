.. _Chap:InputsMultigrid:

Multigrid Inputs
================

Below is a list of the most commonly used multigrid settings options.
To control the nodal projection precede with "nodal_proj", for the MAC projection use "mac_proj", and
for the diffusion solvers use "scalar_diffusion" (tracers, temperature and the component-wise
velocity solve) or "tensor_diffusion" (the tensor velocity solve).  The diffusion solvers read
the verbosity and iteration keys with an "mg_" prefix: mg_verbose, mg_bottom_verbose,
mg_max_iter, mg_bottom_maxiter, mg_maxorder.

+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
|                         |  Description                                                          |   Type      | Default        |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| verbose                 |  Verbosity of multigrid solver                                        |    Int      |   0            |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| bottom_verbose          |  Verbosity of BiCGStab solver                                         |    Int      |   0            |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_rtol                 |  Relative tolerance                                                   |  | Real     |  | 1.e-11      |
|                         |                                                                       |  | float    |  | 1.e-4       |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_atol                 |  Absolute tolerance                                                   |  | Real     |  | 1.e-14      |
|                         |                                                                       |  | float    |  | 1.e-7       |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| maxiter                 |  Maximum number of iterations (diffusion: mg_max_iter)                |    Int      | nodal 100      |
|                         |                                                                       |             | MAC   200      |
|                         |                                                                       |             | diffusion 100  |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| bottom_maxiter          |  Maximum number of iterations in the bottom solver                    |    Int      | nodal 100      |
|                         |  if using bicg, cg, bicgcg or cgbicg (diffusion: mg_bottom_maxiter)   |             | MAC   200      |
|                         |                                                                       |             | diffusion 100  |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_max_coarsening_level |  Maximum number of coarser levels to allow.                           |    Int      |   100          |
|                         |  If set to 0, the bottom solver will be called at the current level   |             |                |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| bottom_solver           |  Which bottom solver to use.                                          |  String     |   bicgcg       |
|                         |  Options are bicgcg, bicgstab, cg, cgbicg, smoother or hypre          |             |                |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+

See AMReX-Hydro's documentation on :ref:`projections inputs <hydro:projections_inputs>` for additional projection options.
