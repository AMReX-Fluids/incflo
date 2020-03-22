Running the MFiX Test Suite
===========================

MFiX-Exa comes with several tests aimed at evaluating software
functionalities. The source files as well as the required input files
for each test are located in the ``tests`` directory. The ``tests``
directory is copied to the build directory during MFiX-Exa configuration
process. When a test is run (see below), output files are stored in
``build_dir/tests/test-name``.

There are various dependencies for comparing test results.

o To compare results to archived flow slices stored in ``AUTOTEST``
directories with the case files, the environment variable ``FEXTRACT``
must point to the location of the AMReX ``fextract`` utility located in
the directory, ``amrex/Tools/PostProcessing/F_Src``. Additionally,
``numdiff`` must be installed for comparing results.

o To compare point-by-point field data, the environment variable
``FCOMPARE`` must point the AMReX utility ``plt_compare_diff_grids``
found in the directory, ``amrex/Tools/PostProcessing/F_Src``.
Additionally, the environment variable ``MFiX_BENCHMARKS_HOME`` must
point to the location of a local regression test data set. See
*Generating local regression test data* for instructions on creating a
local regression test data set.

Run all tests
-------------

.. code:: shell

    > cd to incflo-build-dir
    > ctest

List all tests (without running them)
-------------------------------------

.. code:: shell

    > cd to incflo-build-dir
    > ctest -N

Run a particular test by the index listed in ctest -N
-----------------------------------------------------

.. code:: shell

    > cd to incflo-build-dir
    > ctest -I 3,3             # run the third test

Run a particular test by name
-----------------------------

.. code:: shell

    > cd to incflo-build-dir
    > ctest -R DEM01  # running all tests with "DEM01" in the test name

Run a particular test via make
------------------------------

.. code:: shell

    > cd to incflo-build-dir
    > make run_DEM01-x  # running "DEM01-x" and output to the screen

Run specific
------------

If the environment variable GRID is defined, it specifies which grid
types to run for the test(s). If GRID variable is not defined, the
default is to run the tests for all grid types. > env GRID="tiled" ctest
-R DEM01 # running all tests with "DEM01" for tiled grid > env
GRID="single multiple" ctest -R DEM01 # running all tests with "DEM01"
for single grid and multiple grid > ctest -R DEM01 # running all tests
with "DEM01" for all grid types (single, multiple, tiled)

Run a user-defined case
-----------------------

.. code:: shell

    > ./incflo inputs-myrun

