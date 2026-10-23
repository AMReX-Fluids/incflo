.. _Chap:InputsAlgorithm:

Algorithm Options
=================

The following inputs must be preceded by "incflo."

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
|  advection_type      |  Advection scheme (options: Godunov, MOL, BDS). BDS is not available  |  string     |  Godunov     |
|                      |  with an EB geometry                                                  |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  advect_momentum     | Advect momentum (otherwise use convective differencing)               |  bool       |  false       |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  advect_tracer       |  evolve the tracer equation(s)?                                       |  bool       |  false       |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  trac_is_conservative|  Is tracer conserved? If specified, one entry required per tracer     |  int        |  1           |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  ntrac               |  number of tracers                                                    |  int        |  1           |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  use_temperature     |  include temperature equation? (not available with MOL)               |  bool       |  false       |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  constant_density    |  Only evolve the continuity equation if false                         |  bool       |  true        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  ro_0                |  density (if constant)                                                |  Real       |  1.0         |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  diffusion_type      |  Diffusion type (0 = Explicit, 1 = Crank-Nicolson, 2 = Implicit)      |       int   | 2 (Implicit) |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  mu                  |  viscosity (if constant)                                              |  Real       |  1.0         |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  mu_s                |  scalar diffusivity                                                   |  Real(s)    |  0.0         |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  mu_T                |  thermal conductivity                                                 |  Real       |  0.0         |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  use_divu_source     |  Use the constraint :math:`\nabla \cdot U = S` instead of             |  bool       |  false       |
|                      |  :math:`\nabla \cdot U = 0` (see the Fluid Equations section).        |             |              |
|                      |  Requires advection_type = Godunov and the nodal projection;          |             |              |
|                      |  not yet implemented with embedded boundaries, RZ, direction          |             |              |
|                      |  dependent boundaries, or refinement ratios other than 2 or 4         |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  use_cc_proj         |  Use cell-centered rather than nodal pressure; this changes the       |  bool       | false        |
|                      |  approximate projection we use                                        |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|  use_tensor_solve    |  In velocity solve, use multicomponent :math:`\nabla \cdot \tau`      |  bool       |  true        |
|                      |  otherwise use separate solves for each velocity component            |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

