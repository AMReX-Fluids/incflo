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

The following inputs must be preceded by "bc."

+-----------------+-----------------------------------------------------------------------+-------------+-----------+
|                 | Description                                                           |   Type      | Default   |
+=================+=======================================================================+=============+===========+
| delp_dir        | Direction for specified pressure drop. Note that this direction       |   Int       |   0       |
|                 | should also be periodic.                                              |             |           |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+
| delp            | Pressure drop (Pa)                                                    |   Real      |   0.0     |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+


The following inputs must be preceded by "mfix."

+----------------------+-------------------------------------------------------------------------+----------+-----------+
|                      | Description                                                             |   Type   | Default   |
+======================+=========================================================================+==========+===========+
| geometry             | Which type of EB geometry are we using?                                 |   String |           |
+----------------------+-------------------------------------------------------------------------+----------+-----------+
| levelset__refinement | Refinement factor of levelset resolution relative to level 0 resolution |   Int    | 1         |
+----------------------+-------------------------------------------------------------------------+----------+-----------+
| po_no_par_out        | Let particles exit (default) or bounce-back at pressure outflows        |   Int    | 0         |
+----------------------+-------------------------------------------------------------------------+----------+-----------+
| gravity              | Gravity vector (e.g., mfix.gravity = -9.81  0.0  0.0) [required]        |  Reals   |  None     |
+----------------------+-------------------------------------------------------------------------+----------+-----------+


Setting basic EB walls can be specified by inputs preceded by "xlo", "xhi", "ylo", "yhi", "zlo", and "zhi"

+--------------------+---------------------------------------------------------------------------+-------------+-----------+
|                    | Description                                                               |   Type      | Default   |
+====================+===========================================================================+=============+===========+
| type               | Used to define boundary type. Available options include:                  |  String     |  None     |
|                    |                                                                           |             |           |
|                    | * 'pi'  or 'pressure_inflow'                                              |             |           |
|                    | * 'po'  or 'pressure_outflow'                                             |             |           |
|                    | * 'mi'  or 'mass_inflow'                                                  |             |           |
|                    | * 'nsw' or 'no_slip_wall'                                                 |             |           |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| pressure           | Sets boundary pressure for pressure inflows, outflows and mass inflows    |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| velocity           | Sets boundary velocity for mass inflows                                   |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| location           | Specifies an offset from the domain boundary for no-slip walls            |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+

To specify multiple mass inflows (e.g., define a jet and uniform background flow), provide multiple velocities for the region and define the physical extents of the sub-region. The first velocity is applied to the entire flow plane. Subsequent velocities are successively applied to the specified sub-regions. If multiple sub-regions overlap, the velocity of last specified region is used. An example of a uniform mass inflow with a square-jet centered at (0.5x0.5) is given below.


.. code-block:: none

   xlo.type = "mi"
   xlo.velocity = 0.01  0.10

   xlo.ylo =            0.25
   xlo.yhi =            0.75
   xlo.zlo =            0.25
   xlo.zhi =            0.75


Fluid model settings
--------------------

Enabling the fluid solver and specifying fluid model options.

+----------------------+-------------------------------------------------------------------------+----------+-----------+
|                      | Description                                                             |   Type   | Default   |
+======================+=========================================================================+==========+===========+
| fluid.solve          | Specified name of the fluid or None to disable the fluid solver. The    | String   |  None     |
|                      | name assigned to the fluid solver is used to specify fluid inputs.      |          |           |
+----------------------+-------------------------------------------------------------------------+----------+-----------+


The following inputs must be preceded by the given to the fluid solver e.g., "fluid."

+----------------------+-------------------------------------------------------------------------+----------+-----------+
|                      | Description                                                             |   Type   | Default   |
+======================+=========================================================================+==========+===========+
| density              | Specify which density model to use for fluid [required]                 | String   |  None     |
|                      |   "constant" -- constant density                                        |          |           |
+----------------------+-------------------------------------------------------------------------+----------+-----------+
| density.constant     | Value of constant fluid density [required if density_model= "constant"  |  Real    |  None     |
+----------------------+-------------------------------------------------------------------------+----------+-----------+
| viscosity            | Specify which viscosity model to use for fluid [required]               | String   |  None     |
|                      |   "constant" -- constant viscosity                                      |          |           |
+----------------------+-------------------------------------------------------------------------+----------+-----------+
| viscosity.constant   | Value of constant fluid viscosity [required if viscosity_model="constant"|  Real   |  None     |
+----------------------+-------------------------------------------------------------------------+----------+-----------+

Below is an example for specifying fluid solver model options.

.. code-block:: none

   fluid.solve = myfluid

   myfluid.density = constant
   myfluid.density.constant = 1.0

   myfluid.viscosity = constant
   myfluid.viscosity.constant = 1.8e-5



DEM model settings
------------------

Enabling the DEM solver and specifying model options.

+-------------------------+-------------------------------------------------------------------------+----------+-----------+
|                         | Description                                                             |   Type   | Default   |
+=========================+=========================================================================+==========+===========+
| dem.solve               | Specified name(s) of the DEM types or None to disable the DEM solver.   | String   |  None     |
|                         | The user defined names are used to specify DEM model inputs.            |          |           |
+-------------------------+-------------------------------------------------------------------------+----------+-----------+
| dem.friction_coeff.pp   | Friction coefficient :: particle to particle collisions [required]      | Real     |  None     |
+-------------------------+-------------------------------------------------------------------------+----------+-----------+
| dem.friction_coeff.pw   | Friction coefficient :: particle to wall collisions [required]          | Real     |  None     |
+-------------------------+-------------------------------------------------------------------------+----------+-----------+
| dem.spring_const.pp     | Normal spring constant :: particle to particle collisions [required]    | Real     |  None     |
+-------------------------+-------------------------------------------------------------------------+----------+-----------+
| dem.spring_const.pw     | Normal spring constant :: particle to wall collisions [required]        | Real     |  None     |
+-------------------------+-------------------------------------------------------------------------+----------+-----------+
| dem.spring_tang_fac.pp  | Tangential-to-normal spring constant factor :: particle to particle     | Real     |  None     |
|                         | collisions [required]                                                   |          |           |
+-------------------------+-------------------------------------------------------------------------+----------+-----------+
| dem.spring_tang_fac.pw  | Tangential-to-normal spring constant factor :: particle to wall         | Real     |  None     |
|                         | collisions [required]                                                   |          |           |
+-------------------------+-------------------------------------------------------------------------+----------+-----------+
| dem.damping_tang_fac.pp | Factor relating the tangential damping coefficient to the normal        | Real     |  None     |
|                         | damping coefficient :: particle to particle collisions [required]       |          |           |
+-------------------------+-------------------------------------------------------------------------+----------+-----------+
| dem.damping_tang_fac.pw | Factor relating the tangential damping coefficient to the normal        | Real     |  None     |
|                         | damping coefficient :: particle to wall collisions [required]           |          |           |
+-------------------------+-------------------------------------------------------------------------+----------+-----------+

The following inputs use the DEM type names specified using the `dem.solve` input to define restitution coefficients and
are proceeded with `dem.restitution_coeff`. These must be defined for all solid-solid and solid-wall combinations.

+-------------------------+-------------------------------------------------------------------------+----------+-----------+
|                         | Description                                                             |   Type   | Default   |
+=========================+=========================================================================+==========+===========+
| [solid0].[solid1]       | Specifies the restitution coefficient between solid0 and solid1. Here   | Real     |  None     |
|                         | the order is not important and could be defined as [solid1].[solid0]    |          |           |
+-------------------------+-------------------------------------------------------------------------+----------+-----------+
| [solid0].wall           | Specifies the restitution coefficient between solid0 and the wall.      | Real     |  None     |
|                         | Order is not important and this could be defined as wall.[solid0]       |          |           |
+-------------------------+-------------------------------------------------------------------------+----------+-----------+

Below is an example for specifying the inputs for two DEM solids.

.. code-block:: none

   dem.solve = sand  char

   dem.friction_coeff.pp     =     0.25
   dem.friction_coeff.pw     =     0.15

   dem.spring_const.pp       =   100.0
   dem.spring_const.pw       =   100.0

   dem.spring_tang_fac.pp    =     0.2857
   dem.spring_tang_fac.pw    =     0.2857

   dem.damping_tang_fac.pp   =     0.5
   dem.damping_tang_fac.pw   =     0.5

   dem.restitution_coeff.sand.sand =  0.85
   dem.restitution_coeff.sand.char =  0.88
   dem.restitution_coeff.char.char =  0.90

   dem.restitution_coeff.sand.wall =  0.85
   dem.restitution_coeff.char.wall =  0.89


Region definitions
------------------

Regions are used to define sections of the domain. They may be either boxes, planes or points. They are used in building initial condition regions.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| mfix.regions        | Names given to regions.                                               | String      | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| regions.[region].lo | Low corner of physical region (physical not index space)              |   Reals     | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| regions.[region].hi | High corner of physical region (physical not index space)             |   Reals     | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+

Below is an example for specifying two regions.

.. code-block:: none

   mfix.regions  = full-domain   riser

   regions.full-domain.lo = 0.0000  0.0000  0.0000
   regions.full-domain.hi = 3.7584  0.2784  0.2784

   regions.riser.lo       = 0.0000  0.0000  0.0000
   regions.riser.hi       = 0.1000  0.2784  0.2784



Initial Conditions
------------------

Initial conditions are build from defined regions. The input names are built using the prefix `ic.`, the name of the
region to apply the IC, and the name of the phase (e.g., `myfuild`).

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| ic.regions          | Regions used to define initial conditions.                            | String      | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+

For a fluid phase, the following inputs can be defined.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| volfrac             | Volume fraction [required]                                            | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| pressure            | Fluid pressure                                                        | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| temperature         | Fluid temperature                                                     | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| velocity            | Velocity components                                                   | Reals       | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+


The name of the DEM phases to be defined in the IC region and the packing must be defined.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| ic.[region].solids  | List of solids                                                        | Strings     | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| ic.[region].packing | Specifies how auto-generated particles are placed in the IC region:   | String      | None      |
|                     | * hcp - Hex-centered packing                                          |             |           |
|                     | * random                                                              |             |           |
|                     | * pseudo_random                                                       |             |           |
|                     | * oneper -- one particle per cell                                     |             |           |
|                     | * eightper -- eight particles per cell                                |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+

For each solid, the following inputs may be defined.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| volfrac             | Volume fraction                                                       | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| temperature         | Fluid temperature                                                     | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| velocity            | Velocity components                                                   | Reals       | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| diameter            | Method to specify particle diameter in the IC region. This is         | String      | None      |
|                     | only used for auto-generated particles.                               |             |           |
|                     | * constant  -- specified constant                                     |             |           |
|                     | * uniform   -- uniform distribution                                   |             |           |
|                     | * normal    -- normal distribution                                    |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| diameter.constant   | Value of specified constant particle density                          | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| diameter.mean       | Distribution mean                                                     | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| diameter.std        | Distribution standard deviation                                       | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| diameter.min        | Minimum diameter to clip distribution                                 | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| diameter.max        | Maximum diameter to clip distribution                                 | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| density             | Method to specify particle density in the IC region. This is          | String      | None      |
|                     | only used for auto-generated particles.                               |             |           |
|                     | * constant  -- specified constant                                     |             |           |
|                     | * uniform   -- uniform distribution                                   |             |           |
|                     | * normal    -- normal distribution                                    |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| density.constant    | Value of specified constant particle density                          | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| density.mean        | Distribution mean                                                     | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| density.std         | Distribution standard deviation                                       | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| density.min         | Minimum density to clip distribution                                  | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| density.max         | Maximum density to clip distribution                                  | Real        | None      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+


Below is an example for specifying an initial condition for a fluid (fluid) and one DEM solid (solid0).

.. code-block:: none

   ic.regions  = bed

   ic.bed.fluid.volfrac   =  0.725

   ic.bed.fluid.velocity  =  0.015  0.00  0.00

   ic.bed.solids  = solid0
   ic.bed.packing = pseudo_random

   ic.bed.solid0.volfrac  =  0.275

   ic.bed.solid0.velocity =  0.00  0.00  0.00

   ic.bed.solid0.diameter = constant
   ic.bed.solid0.diameter.constant =  100.0e-6

   ic.bed.solid0.density  = constant
   ic.bed.solid0.density.constant  = 1000.0
