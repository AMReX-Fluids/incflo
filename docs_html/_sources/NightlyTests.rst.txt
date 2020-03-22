.. _Chap:NightlyTesting :

Nightly Tests
=============

The following regression tests are run nightly with MFiX-Exa.   The plotfiles generated in each night's test 
are compared with the benchmark plotfiles using the AMReX :cpp:`fcompare` utility to compare the mesh data
and :cpp:`particle_compare` to compare the particle data.

The results of these tests can be found at https://ccse.lbl.gov/pub/RegressionTesting/MFIX-Exa/

Below Ng = number of grids, Npa = number of particles, Np = number of MPI ranks.

"Auto" means the particles were generated automatically with the random number
generator; if "Auto" is not specified the particle data were read in from "particle_input.dat"

These first tests have both fluid and particles and are run in rectangular geometries;
all tests except DEM06 use drag type "BVK2". 

"NSW" means "No Slip Wall" and "Per" is "periodic."
"MI/PO" refers to Mass Inflow at the low end of the domain and Pressure Outflow at the high end.

+-------------------+----+--------+------+-------+----+----+----------------------+
| Test              | nx | bc_x   | EB   | Npa   | Ng | Np | What does this test? |
|                   | ny | bc_y   |      |       |    |    |                      |
|                   | nz | bc_z   |      |       |    |    |                      |
+===================+====+========+======+=======+====+====+======================+
| BENCH01           | 32 | Per    | None | 5005  | 1  | 1  | Triply periodic      |
| Size0001          | 32 | Per    |      |       |    |    |                      |
|                   | 32 | Per    |      |       |    |    |                      |
+-------------------+----+--------+------+-------+----+----+----------------------+
| BENCH01           | 64 | Per    | None | 40040 | 8  | 4  | Replicate            |
| Size0001          | 64 | Per    |      |       |    |    |                      |
| replicate         | 64 | Per    |      |       |    |    |                      |
+-------------------+----+--------+------+-------+----+----+----------------------+
| BENCH01           | 32 | Per    | None | 5005  | 8  | 4  | Restart              |
| Size0001          | 32 | Per    |      |       |    |    |                      |
| restart           | 32 | Per    |      |       |    |    |                      |
+-------------------+----+--------+------+-------+----+----+----------------------+
| BENCH02           | 10 | Per    | None | 1611  | 1  | 1  | Mixed NSW / Per      |
| Size0001          | 10 | NSW    |      |       |    |    |                      |
|                   | 10 | Per    |      |       |    |    |                      |
+-------------------+----+--------+------+-------+----+----+----------------------+
| BENCH02           | 10 | NSW    | None | 1611  | 1  | 1  | NSW on all faces     |
| Size0001          | 10 | NSW    |      |       |    |    |                      |
| walls             | 10 | NSW    |      |       |    |    |                      |
+-------------------+----+--------+------+-------+----+----+----------------------+
| BENCH03           | 4  | Per    | None | 2500  | 1  | 1  | Mixed MI/PO + Per    |
| Size0001          | 50 | MI/PO  |      |       |    |    |                      |
|                   | 4  | Per    |      |       |    |    |                      |
+-------------------+----+--------+------+-------+----+----+----------------------+
| BENCH04           | 4  | Per    | None | 224   | 1  | 1  | Triply periodic      |
| Size0001          | 50 | Per    |      |       |    |    |                      |
|                   | 4  | Per    |      |       |    |    |                      |
+-------------------+----+--------+------+-------+----+----+----------------------+
| DEM06             | 5  | Per    | None | 1     | 10 | 4  | Single particle      |
| z multiple        | 5  | Per    |      |       |    |    | falling in fluid     |
|                   | 50 | MI/PO  |      |       |    |    | (user_drag)          |
+-------------------+----+--------+------+-------+----+----+----------------------+

This second set of tests have both fluid and particles and are run in cylindrial geometries
interior to the domain boundaries; they also use drag type "BVK2".  Here "IGN" means
those domain boundaries should be ignored because they are outside the EB boundary.

+-------------------+----+-------+------+--------+----+----+----------------------+
| Test              | nx | bc_x  | EB   | Npa    | Ng | Np | What does this test? |
|                   | ny | bc_y  |      |        |    |    |                      |
|                   | nz | bc_z  |      |        |    |    |                      |
+===================+====+=======+======+========+====+====+======================+
| BENCH05           | 40 | MI/PO | Cyl  | 7949   | 4  | 4  | EB in parallel       |
| Size0008          | 10 | IGN   |      | Auto   |    |    |                      |
|                   | 10 | IGN   |      |        |    |    |                      |
+-------------------+----+-------+------+--------+----+----+----------------------+
| BENCH05           | 40 | MI/PO | Cyl  | 7968   | 4  | 1  | EB in serial         |
| Size0008          | 10 | IGN   |      | Auto   |    |    |                      |
| serial            | 10 | IGN   |      |        |    |    |                      |
+-------------------+----+-------+------+--------+----+----+----------------------+
| BENCH05           | 40 | MI/PO | Cyl  | 36672  | 16 | 4  | Regrid & dual grid   |
| Size0008          | 20 | IGN   |      | Auto   |    |    |                      |
| medium            | 20 | IGN   |      |        |    |    |                      |
+-------------------+----+-------+------+--------+----+----+----------------------+
| BENCH05           | 40 | MI/PO | Cyl  | 157106 | 16 | 4  | Regrid & dual grid   |
| Size0008          | 40 | IGN   |      | Auto   |    |    |                      |
| wide              | 40 | IGN   |      |        |    |    |                      |
+-------------------+----+-------+------+--------+----+----+----------------------+
| BENCH06           | 40 | Per   | Cyl  | 627    | 4  | 1  | EB                   |
| Size0008          | 10 | IGN   |      | Auto   |    |    | with periodic        |
| serial            | 10 | IGN   |      |        |    |    | serial               |
+-------------------+----+-------+------+--------+----+----+----------------------+
| BENCH06           | 40 | Per   | Cyl  | 624    | 4  | 4  | EB                   |
| Size0008          | 10 | IGN   |      | Auto   |    |    | with periodic        |
|                   | 10 | IGN   |      |        |    |    | parallel             |
+-------------------+----+-------+------+--------+----+----+----------------------+


This third set of tests is particles-only in rectangular geometries.

+-------------------+----+-------+------+--------+----+----+----------------------+
| Test              | nx | bc_x  | EB   | Npa    | Ng | Np | What does this test? |
|                   | ny | bc_y  |      |        |    |    |                      |
|                   | nz | bc_z  |      |        |    |    |                      |
+===================+====+=======+======+========+====+====+======================+
| DEM01             | 4  | NSW   | None |   1    | 1  | 1  | Particle  only       |
| x single          | 4  | NSW   |      |        |    |    |                      |
|                   | 4  | NSW   |      |        |    |    |                      |
+-------------------+----+-------+------+--------+----+----+----------------------+
| DEM03             | 5  | Per   | None |   2    | 1  | 1  | Particles only       |
| z single          | 5  | Per   |      |        |    |    |                      |
|                   | 2  | NSW   |      |        |    |    |                      |
+-------------------+----+-------+------+--------+----+----+----------------------+
| DEM04             | 4  | NSW   | None |   1    | 1  | 1  | Particles only       |
| z single          | 4  | Per   |      |        |    |    |                      |
|                   | 4  | Per   |      |        |    |    |                      |
+-------------------+----+-------+------+--------+----+----+----------------------+
