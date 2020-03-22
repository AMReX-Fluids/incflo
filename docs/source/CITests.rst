.. _Chap:CITesting :

Continuous Integration
======================

The following regression tests are run every time a commit is pushed to the main
MFiX-Exa repository on the NETL gitlab.

For each of the tests in the chart below, there are
three directional variations; these are identified in the repository as, 
for example, FLD01-x, FLD01-y, and FLD01-z.  

For each direction, where appropriate, there are multiple versions, with the following notations:

  * SGS: single grid serial

  * MGS: multiple grid serial

  * TGS: tiled grid serial

  * MGP: multiple grid parallel

Below Ng = number of grids, Npa = number of particles, Np = number of MPI ranks.

All the FLD cases are fluid-only and steady state.
All the DEM cases are particle-only except for DEM06 and DEM07 which are fluid and particles;
these both use the "BVK2" drag type.
In all cases the particle data were read in from "particle_input.dat"

None of these tests have non-rectangular geometry.

"NSW" means "No Slip Wall" and "Per" is "periodic."
"MI/PO" refers to Mass     Inflow at the low end of the domain and Pressure Outflow at the high end.
"PI/PO" refers to Pressure Inflow at the low end of the domain and Pressure Outflow at the high end.

Additional detail about these problems is given in tests/README.md

Single-grid, single-process (SGS) particle-only tests:

+-------+----+----+----+------+------+-------+-----+--------------------+
| Test  | nx | ny | nz | bc_x | bc_y | bc_z  | Npa | Description        |
+=======+====+====+====+======+======+=======+=====+====================+
| DEM01 | 2  | 5  | 5  | NSW  | Per  | Per   | 1   | Freely falling     |
|       |    |    |    |      |      |       |     | particle with      |
|       |    |    |    |      |      |       |     | wall collision     |
+-------+----+----+----+------+------+-------+-----+--------------------+
| DEM02 | 2  | 5  | 5  | NSW  | Per  | Per   | 1   | Multiple bounces   |
|       |    |    |    |      |      |       |     | with bounce height |
|       |    |    |    |      |      |       |     | measured           |
+-------+----+----+----+------+------+-------+-----+--------------------+
| DEM03 | 2  | 5  | 5  | NSW  | Per  | Per   | 2   | Two stacked        |
|       |    |    |    |      |      |       |     | compressed         |
|       |    |    |    |      |      |       |     | particles          |
+-------+----+----+----+------+------+-------+-----+--------------------+
| DEM04 | 4  | 4  | 4  | NSW  | Per  | Per   | 1   | Single particle    |
|       |    |    |    |      |      |       |     | slipping on a      |
|       |    |    |    |      |      |       |     | rough surface      |
+-------+----+----+----+------+------+-------+-----+--------------------+
| DEM05 | 5  |  2 |  5 | Per  | Per  | Per   | 93  | Oblique particle   |
|       |    |    |    |      |      |       |     | collisions         |
|       |    |    |    |      |      |       |     |                    |
+-------+----+----+----+------+------+-------+-----+--------------------+


Steady-state fluid-only tests:

+-------+-----+----+----+----+-------+------+------+-----+----------------------+
| Test  |     | nx | ny | nz | bc_x  | bc_y | bc_z | Ng  | Np                   |
+=======+=====+====+====+====+=======+======+======+=====+======================+
| FLD01 |     | 8  | 8  | 4  | Per   | NSW  | Per  | Poiseuille flow            |
+-------+-----+----+----+----+-------+------+------+-----+----+-----------------+
|       | SGS |    |    |    |       |      |      | 1   | 1  |                 |
|       | MGS |    |    |    |       |      |      | 4   | 1  |                 |
|       | MGP |    |    |    |       |      |      | 4   | 8  |                 |
+-------+-----+----+----+----+-------+------+------+-----+----+-----------------+
| FLD02 |     | 80 |16  | 16 | MI/PO | NSW  | NSW  | Couette flow               |
+-------+-----+----+----+----+-------+------+------+-----+----+-----------------+
|       | SGS |    |    |    |       |      |      | 1   | 1  |                 |
|       | MGS |    |    |    |       |      |      | 40  | 1  |                 | 
|       | MGP |    |    |    |       |      |      | 40  | 8  |                 |
+-------+-----+----+----+----+-------+------+------+-----+----+-----------------+
| FLD03 |     | 8  | 8  | 4  | PI/PO | NSW  | Per  | Poiseuille flow            |
+-------+-----+----+----+----+-------+------+------+-----+----+-----------------+
|       | SGS |    |    |    |       |      |      | 1   | 1  |                 |
|       | MGS |    |    |    |       |      |      | 4   | 1  |                 |
|       | MGP |    |    |    |       |      |      | 4   | 8  |                 |
+-------+-----+----+----+----+-------+------+------+-----+----+-----------------+

Coupled particle/fluid tests:

+-------+-----+----+----+----+------+------+------+------+----+--------------------+
| Test  |     | nx | ny | nz | bc_x | bc_y | bc_z | Npa  | Ng | Np                 |
+=======+=====+====+====+====+======+======+======+======+====+====================+
| DEM06 |     | 50 | 5  | 5  |  NSW | NSW  | NSW  | 1    | Single particle falling | 
|       |     |    |    |    |      |      |      |      | under gravity           |
+-------+-----+----+----+----+------+------+------+------+----+----+---------------+
|       | SGS |    |    |    |      |      |      |      | 1  | 1  |               |
|       | MGS |    |    |    |      |      |      |      | 10 | 1  |               |
|       | MGP |    |    |    |      |      |      |      | 10 | 8  |               |
+-------+-----+----+----+----+------+------+------+------+----+----+---------------+
| DEM07 |     | 20 | 20 | 20 | Per  | Per  | Per  | 1222 | Homogeneous cooling     |
|       |     |    |    |    |      |      |      |      | system                  |
+-------+-----+----+----+----+------+------+------+------+----+----+---------------+
|       | SGS |    |    |    |      |      |      |      | 1  | 1  |               |
|       | MGS |    |    |    |      |      |      |      | 8  | 1  |               |
|       | MGP |    |    |    |      |      |      |      | 8  | 8  |               |
+-------+-----+----+----+----+------+------+------+------+----+----+---------------+
