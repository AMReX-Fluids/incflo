.. _Chap:InputsLoadBalancing:

Gridding and Load Balancing
===========================

The following inputs must be preceded by "amr" and determine how we create the grids and how often we regrid.

+----------------------+-----------------------------------------------------------------------+-------------+-----------+
|                      | Description                                                           |   Type      | Default   |
+======================+=======================================================================+=============+===========+
| regrid_int           | How often to regrid (in number of steps at level 0)                   |   Int       |    -1     |
|                      | if regrid_int = -1 then no regridding will occur                      |             |           |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_x      | Maximum number of cells at level 0 in each grid in x-direction        |    Int      | 32        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_y      | Maximum number of cells at level 0 in each grid in y-direction        |    Int      | 32        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_z      | Maximum number of cells at level 0 in each grid in z-direction        |    Int      | 32        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_x    | Each grid must be divisible by blocking_factor_x in x-direction       |    Int      |  8        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_y    | Each grid must be divisible by blocking_factor_y in y-direction       |    Int      |  8        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_z    | Each grid must be divisible by blocking_factor_z in z-direction       |    Int      |  8        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+

The following inputs must be preceded by "fabarray_mfiter" and determine how we create the logical tiles:

+----------------------+-----------------------------------------------------------------------+----------+-------------+
|                      | Description                                                           | Type     | Default     |
+======================+=======================================================================+==========+=============+
| tile_size            | Maximum number of cells in each direction for (logical) tiles         | IntVect  | 1024000     |
|                      |        (3D CPU-only)                                                  |          | 1024000,8,8 |
+----------------------+-----------------------------------------------------------------------+----------+-------------+

The following inputs must be preceded by "particles"

+-------------------+-----------------------------------------------------------------------+-------------+--------------+
|                   | Description                                                           |   Type      | Default      |
+===================+=======================================================================+=============+==============+
| max_grid_size_x   | Maximum number of cells at level 0 in each grid in x-direction        |    Int      | 32           |
|                   | for grids in the ParticleBoxArray if dual_grid is true                |             |              |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| max_grid_size_y   | Maximum number of cells at level 0 in each grid in y-direction        |    Int      | 32           |
|                   | for grids in the ParticleBoxArray if dual_grid is true                |             |              |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| max_grid_size_z   | Maximum number of cells at level 0 in each grid in z-direction        |    Int      | 32           |
|                   | for grids in the ParticleBoxArray if dual_grid is true.               |             |              |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| tile_size         | Maximum number of cells in each direction for (logical) tiles         |  IntVect    | 1024000,8,8  |
|                   | in the ParticleBoxArray if dual_grid is true.                         |             |              |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+

Note that when running a granular simulation, i.e., no fluid phase, :cpp:`mfix.dual_grid` must be 0. Hence, 
the :cpp:`particles.max_grid_size` (in each direction) have no meaning. Therefore the fluid grid and tile 
sizes should be set for particle load balancing. It may also be necessary to set the blocking factors to 1. 


The following inputs must be preceded by "mfix" and determine how we load balance:

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| dual_grid            | If true then use the "dual_grid" approach for load balancing          |  Bool       | False        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| load_balance_fluid   | Only relevant if (dual_grid); if so do we also regrid mesh data       |  Int        | 1            |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| load_balance_type    | What strategy to use for load balancing                               |  String     | KnapSack     |
|                      | Options are "KnapSack"or "SFC"                                        |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| knapsack_weight_type | What weighting function to use if using Knapsack load balancing       |  String     | RunTimeCosts |
|                      | Options are "RunTimeCosts" or "NumParticles""                         |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| knapsack_nmax        | Maximum number of grids per MPI process if using knapsack algorithm   |  Int        | 128          | 
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

