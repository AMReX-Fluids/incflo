.. _Chap:InputsLoadBalancing:

Gridding and Load Balancing Inputs
==================================

The following inputs must be preceded by "amr" and determine how we create the grids and how often we regrid. The most
commonly used options are listed here. See the AMReX documentation on :ref:`amrex:chap:inputs` for additional
options.

+----------------------+-----------------------------------------------------------------------+-------------+---------------+
|                      | Description                                                           |   Type      | Default       |
+======================+=======================================================================+=============+===============+
| regrid_int           | How often to regrid (in number of steps at level 0)                   |   Int       |    -1         |
|                      | if regrid_int = -1 then no regridding will occur                      |             |               |
+----------------------+-----------------------------------------------------------------------+-------------+---------------+
| max_grid_size_x      | Maximum number of cells in each grid in x-direction (on all levels,   |    Int      | | 2D: 128     |
|                      | unless one value per level is given)                                  |             | | 3D CPU: 32  |
|                      |                                                                       |             | | 3D GPU: 64  |
+----------------------+-----------------------------------------------------------------------+-------------+---------------+
| max_grid_size_y      | Maximum number of cells in each grid in y-direction (on all levels,   |    Int      | | 2D: 128     |
|                      | unless one value per level is given)                                  |             | | 3D CPU: 32  |
|                      |                                                                       |             | | 3D GPU: 64  |
+----------------------+-----------------------------------------------------------------------+-------------+---------------+
| max_grid_size_z      | Maximum number of cells in each grid in z-direction (3D) (on all      |    Int      | | CPU: 32     |
|                      | levels, unless one value per level is given)                          |             | | GPU: 64     |
+----------------------+-----------------------------------------------------------------------+-------------+---------------+
| blocking_factor_x    | Each grid must be divisible by blocking_factor_x in x-direction       |    Int      |  8            |
+----------------------+-----------------------------------------------------------------------+-------------+---------------+
| blocking_factor_y    | Each grid must be divisible by blocking_factor_y in y-direction       |    Int      |  8            |
+----------------------+-----------------------------------------------------------------------+-------------+---------------+
| blocking_factor_z    | Each grid must be divisible by blocking_factor_z in z-direction (3D)  |    Int      |  8            |
+----------------------+-----------------------------------------------------------------------+-------------+---------------+

The following inputs must be preceded by "fabarray" and determine how we create the logical tiles:

+----------------------+-----------------------------------------------------------------------+----------+-------------+
|                      | Description                                                           | Type     | Default     |
+======================+=======================================================================+==========+=============+
| mfiter_tile_size     | Maximum number of cells in each direction for (logical) tiles         | IntVect  | 1024000     |
|                      |                                                                       |          |             |
|                      | (3D CPU-only)                                                         |          | 1024000,8,8 |
+----------------------+-----------------------------------------------------------------------+----------+-------------+
