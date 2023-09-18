.. _Chap:InputsMultigrid:

Multigrid Inputs
================

Below is a list of the most commonly used multigrid settings options.
To control the nodal projection precede with "nodal_proj", for the MAC projection use "mac_proj", and
for the diffusion solver use "diffusion"

+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
|                         |  Description                                                          |   Type      | Default        |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| verbose                 |  Verbosity of multigrid solver                                        |    Int      |   0            |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| bottom_verbose          |  Verbosity of BiCGStab solver                                         |    Int      |   0            |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| rtol                    |  Relative tolerance                                                   |  | Real     |  | 1.e-11      |
|                         |                                                                       |  | float    |  | 1.e-4       |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| atol                    |  Absolute tolerance                                                   |  | Real     |  | 1.e-14      |
|                         |                                                                       |  | float    |  | 1.e-7       |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| maxiter                 |  Maximum number of iterations                                         |    Int      | nodal 100      |
|                         |                                                                       |             | MAC   200      |
|                         |                                                                       |             | diffusion 100  |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| bottom_maxiter          |  Maximum number of iterations in the                                  |    Int      | nodal 100      |
|                         |  bottom solver if using bicg, cg, bicgcg or cgbicg                    |             | MAC   200      |
|                         |                                                                       |             | diffusion 100  |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_max_coarsening_level |  Maximum number of coarser levels to allow.                           |    Int      |   100          |
|                         |  If set to 0, the bottom solver will be called at the current level   |             |                |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| bottom_solver           |  Which bottom solver to use.                                          |  String     |   bicgcg       |
|                         |  Options are bicgcg, bicgstab, cg, cgbicg, smoother or hypre          |             |                |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+

See AMReX-Hydro's documentation on :ref:`projections inputs <hydro:projections_inputs>` for additional projection options.
