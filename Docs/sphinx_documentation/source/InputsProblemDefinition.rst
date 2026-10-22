Problem Definition
==================

The following inputs must be preceded by "amr."

+-------------------+-----------------------------------------------------------------------+-------------+-----------+
|                   | Description                                                           |   Type      | Default   |
+===================+=======================================================================+=============+===========+
| n_cell            | Number of cells at level 0 in each coordinate direction               | Int Int Int | None      |
+-------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_level         | Maximum level of refinement allowed (0 when single-level)             |    Int      | None      |
+-------------------+-----------------------------------------------------------------------+-------------+-----------+

The following inputs must be preceded by "geometry."

+-----------------+-----------------------------------------------------------------------+-------------+-----------+
|                 | Description                                                           |   Type      | Default   |
+=================+=======================================================================+=============+===========+
| coord_sys       | 0 for Cartesian                                                       |   Int       |   0       |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+
| is_periodic     | 1 for true, 0 for false (one value for each coordinate direction)     |   Ints      | 0 0 0     |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+
| prob_lo         | Low corner of physical domain (physical not index space)              |   Reals     | None      |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+
| prob_hi         | High corner of physical domain (physical not index space)             |   Reals     | None      |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+


The following inputs must be preceded by "incflo."

+----------------------+-------------------------------------------------------------------------+----------+-----------+
|                      | Description                                                             |   Type   | Default   |
+======================+=========================================================================+==========+===========+
| geometry             | Which type of EB geometry are we using? (EB builds only)                |   String | None      |
|                      | Options are box, cylinder, annulus, sphere, jcap, stl, chkptfile,       |          |           |
|                      | csg (CSG builds), and, in 3D only, twocylinders, spherecube and         |          |           |
|                      | tuscan. If not set (or set to all_regular), there is no EB.             |          |           |
+----------------------+-------------------------------------------------------------------------+----------+-----------+
| gravity              | Gravity vector (e.g., incflo.gravity = -9.81  0.0  0.0)                 |  Reals   | (0, 0, 0) |
+----------------------+-------------------------------------------------------------------------+----------+-----------+
| delp                 | Pressure drop across the domain in each coordinate direction, which     |  Reals   | (0, 0, 0) |
|                      | drives the flow (see :ref:`sec:BackgroundPressure`). At most one value  |          |           |
|                      | can be non-zero, and not in a direction with gravity or with a pressure |          |           |
|                      | inflow/outflow pair                                                     |          |           |
+----------------------+-------------------------------------------------------------------------+----------+-----------+


Setting basic boundary conditions can be specified by inputs preceded by "xlo", "xhi", "ylo", "yhi", "zlo", and "zhi".
No boundary type is given for the faces of a periodic direction (``geometry.is_periodic``).

+--------------------+---------------------------------------------------------------------------+-------------+-----------+
|                    | Description                                                               |   Type      | Default   |
+====================+===========================================================================+=============+===========+
| type               | Used to define boundary type. Available options include:                  |  String     |  None     |
|                    |                                                                           |             |           |
|                    | * 'pi'  or 'pressure_inflow'                                              |             |           |
|                    | * 'po'  or 'pressure_outflow'                                             |             |           |
|                    | * 'mi'  or 'mass_inflow'                                                  |             |           |
|                    | * 'dd'  or 'direction_dependent'                                          |             |           |
|                    | * 'nsw' or 'no_slip_wall'                                                 |             |           |
|                    | * 'sw' or 'slip_wall'                                                     |             |           |
|                    | * 'mixed'                                                                 |             |           |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| pressure           | Sets boundary pressure for pressure inflows, pressure outflows and mixed  |    Real     |  None     |
|                    | boundaries (required for these types)                                     |             |           |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| velocity           | Sets Dirichlet boundary condition for velocity at mass inflow,            |    Reals    | (0, 0, 0) |
|                    | direction-dependent and mixed boundaries. At a no-slip wall, the          |             |           |
|                    | tangential components give the velocity of the wall                       |             |           |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| density            | Sets Dirichlet boundary condition for density at mass inflow,             |    Real     |  1.0      |
|                    | direction-dependent and mixed boundaries                                  |             |           |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| tracer             | Sets Dirichlet boundary condition for tracer (one value per tracer) at    |    Reals    |  0.0      |
|                    | mass inflow, direction-dependent and mixed boundaries. At a wall, the     |             |           |
|                    | tracer is Dirichlet if this is given, and has zero flux otherwise         |             |           |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| temperature        | Sets Dirichlet boundary condition for temperature, as for tracer          |    Real     |  1.0      |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+

   The 'mixed' boundary type allows for inflow and outflow on the same domain face.
   The implementation requires that there is EB separating the inflow and outflow regions,
   requires ``incflo.advection_type = Godunov``,
   and only currently treats constant viscosity (i.e. requires the inputs file contains  ``incflo.use_tensor_solve = false``).
   It does allow for non-constant scalar diffusivity. To create a new problem setup with mixed BCs, one must additionally
   specify the Dirichlet and Neumann areas in ``prob/prob_bc.cpp``. In this file there are functions to set the needed BC info for the diffusion solver, the MAC projection, and with the same function, the nodal projection and advection. Comments in the code provide a detailed explanation for each instance. For additional details on mixed BCs, also see AMReX-Hydro's documentation (:ref:`hydro:bcs`).

..
   To create a new setup:
   initial conditions go in prob/prob_init_fluid.cpp
   inflow boundary conditions are set in prob/prob_bc.H -- the BCType for mixed bc is set at foextrap here to allow the outflow to get appropriately filled, which would already have happened before reaching these functions.
