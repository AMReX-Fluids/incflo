.. _Chap:InputsMultigrid:

Multigrid Inputs
================

incflo uses AMReX's multigrid solver, MLMG, for the projections and for the implicit diffusion
solves. The AMReX documentation on :ref:`amrex:chap:linearsolvers` describes the solver and its
parameters (e.g. the available bottom solvers).

Projections
-----------

The MAC and nodal projections are done by AMReX-Hydro's ``MacProjector`` and ``NodalProjector``
classes, which read their own options (verbosity, maximum number of iterations, bottom solver, ...)
with the prefix "mac_proj" or "nodal_proj". These are documented in AMReX-Hydro's documentation on
:ref:`projections inputs <hydro:projections_inputs>`. The cell-centered projection
(``incflo.use_cc_proj = 1``) is also done with the ``MacProjector``, so it uses the "mac_proj"
options.

In addition, incflo reads the following, which must be preceded by "mac_proj" or "nodal_proj":

+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
|                         |  Description                                                          |   Type      | Default        |
+=========================+=======================================================================+=============+================+
| mg_rtol                 |  Relative tolerance                                                   |  | Real     |  | 1.e-11      |
|                         |                                                                       |  | float    |  | 1.e-4       |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_atol                 |  Absolute tolerance                                                   |  | Real     |  | 1.e-14      |
|                         |                                                                       |  | float    |  | 1.e-7       |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_max_coarsening_level |  Maximum number of coarser levels to allow.                           |    Int      |   100          |
|                         |  If set to 0, the bottom solver will be called at the current level   |             |                |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+

Diffusion
---------

The diffusion solves are set up by incflo itself. Their options must be preceded by
"scalar_diffusion" (tracers, temperature and, with ``incflo.use_tensor_solve = false``, the
component-wise velocity solve) or "tensor_diffusion" (the velocity solve with
``incflo.use_tensor_solve = true``):

+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
|                         |  Description                                                          |   Type      | Default        |
+=========================+=======================================================================+=============+================+
| verbose                 |  Verbosity of incflo's diffusion routines                             |    Int      |   0            |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_verbose              |  Verbosity of the multigrid solver                                    |    Int      |   0            |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_bottom_verbose       |  Verbosity of the bottom solver ("scalar_diffusion" also accepts      |    Int      |   0            |
|                         |  bottom_verbose)                                                      |             |                |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_max_iter             |  Maximum number of multigrid iterations                               |    Int      |   100          |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_bottom_maxiter       |  Maximum number of iterations in the bottom solver                    |    Int      |   100          |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_max_fmg_iter         |  Number of full multigrid cycles before switching to V-cycles         |    Int      |   0            |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_max_coarsening_level |  Maximum number of coarser levels to allow.                           |    Int      |   100          |
|                         |  If set to 0, the bottom solver will be called at the current level   |             |                |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_maxorder             |  Maximum order (2, 3 or 4) of the interpolation used at boundaries    |    Int      |   2            |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_rtol                 |  Relative tolerance                                                   |  | Real     |  | 1.e-11      |
|                         |                                                                       |  | float    |  | 1.e-4       |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| mg_atol                 |  Absolute tolerance                                                   |  | Real     |  | 1.e-14      |
|                         |                                                                       |  | float    |  | 1.e-7       |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| bottom_solver           |  Which bottom solver to use: smoother or hypre. Any other value       |  String     |   bicgstab     |
|                         |  gives MLMG's default bottom solver for these operators, BiCGStab     |             |                |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| num_pre_smooth          |  Number of smoother iterations when going down the V-cycle            |    Int      |   2            |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| num_post_smooth         |  Number of smoother iterations when going up the V-cycle              |    Int      |   2            |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+
| use_gauss_seidel        |  "tensor_diffusion" only: use Gauss-Seidel rather than Jacobi         |    Bool     |   true         |
|                         |  smoothing (not used when there is an EB)                             |             |                |
+-------------------------+-----------------------------------------------------------------------+-------------+----------------+

The hypre bottom solver requires incflo to be built with HYPRE support (``USE_HYPRE = TRUE`` or
``-DINCFLO_HYPRE=yes``). hypre's own options take the prefix "hypre"; see the
`External Solvers <https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html#external-solvers>`_
section of the AMReX documentation.
