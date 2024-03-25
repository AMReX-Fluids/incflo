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
| geometry             | Which type of EB geometry are we using?                                 |   String |           |
+----------------------+-------------------------------------------------------------------------+----------+-----------+
| gravity              | Gravity vector (e.g., incflo.gravity = -9.81  0.0  0.0)                 |  Reals   | (0, 0, 0) |
+----------------------+-------------------------------------------------------------------------+----------+-----------+
| delp                 | Pressure drop (Pa)                                                      |   Real   | (0, 0, 0) |
+----------------------+-------------------------------------------------------------------------+----------+-----------+


Setting basic boundary conditions can be specified by inputs preceded by "xlo", "xhi", "ylo", "yhi", "zlo", and "zhi"

+--------------------+---------------------------------------------------------------------------+-------------+-----------+
|                    | Description                                                               |   Type      | Default   |
+====================+===========================================================================+=============+===========+
| type               | Used to define boundary type. Available options include:                  |  String     |  None     |
|                    |                                                                           |             |           |
|                    | * 'pi'  or 'pressure_inflow'                                              |             |           |
|                    | * 'po'  or 'pressure_outflow'                                             |             |           |
|                    | * 'mi'  or 'mass_inflow'                                                  |             |           |
|                    | * 'nsw' or 'no_slip_wall'                                                 |             |           |
|                    | * 'mixed'                                                                 |             |           |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| pressure           | Sets boundary pressure for pressure inflows, outflows and mass inflows    |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| velocity           | Sets boundary velocity for mass inflows                                   |    Reals    | (0, 0, 0) |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| density            | Sets boundary density for mass inflows                                    |    Real     |  1.0      |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| tracer             | Sets boundary tracer for mass inflows                                     |    Real     |  0.0      |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+

   The 'mixed' boundary type allows for inflow and outflow on the same domain face.
   The implementation requires that there is EB separating the inflow and outflow regions,
   and only currently treats constant viscosity (i.e. requires the inputs file contains  ``incflo.use_tensor_solve = false``).
   It does allow for non-constant scalar diffusivity. To create a new problem setup with mixed BCs, one must additionally
   specify the Dirchlet and Neumann areas in ``prob/prob_bc.cpp`` In this file there are functions to set the needed BC info for the diffusion solver, the MAC projection, and with the same function, the nodal projection and advection. Comments in the code provide a detailed explaination for each instance. For addditional details on mixed BCs, also see AMReX-Hydro's documentation (:ref:`bcs`).

..
   To create a new setup:
   initial conditions go in prob/prob_init_fluid.cpp
   inflow boundary conditions are set in prob/prob_bc.H -- the BCType for mixed bc is set at foextrap here to allow the outflow to get appropriately filled, which would already have happened before reaching these functions.

