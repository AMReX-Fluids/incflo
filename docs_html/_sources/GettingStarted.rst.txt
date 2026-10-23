.. _Chap:GettingStarted:

Getting Started
***************
This section walks you through a brief introduction to using incflo.

..
    alternative to `` `` is :code:` `


.. _sec:Download:

Downloading the code
====================
incflo is built on top of the AMReX framework.  In order to run
incflo, you must download separate git modules for incflo, AMReX
and AMReX-Hydro.

First, make sure that git is installed on your machine.

#. Download the AMReX repository by typing:

   .. code:: shell

             git clone https://github.com/AMReX-Codes/amrex.git

   This will create a folder called ``amrex/`` on your machine.
   For additional information on downloading the AMReX source code see `AMReX user guide - Downloading the code <https://amrex-codes.github.io/amrex/docs_html/GettingStarted.html>`_.


#. Download the AMReX-Hydro repository by typing:

   .. code:: shell

         git clone https://github.com/AMReX-Fluids/AMReX-Hydro.git

   This will create a folder called ``AMReX-Hydro/`` on your machine.

#. Download the incflo repository by typing:

   .. code:: shell

         git clone https://github.com/AMReX-Fluids/incflo.git

   This will create a folder called ``incflo/`` on your machine.

You will want to periodically update each of these repositories
by typing ``git pull`` within each repository.

Next, you can choose to either :ref:`sec:GNUMake` or :ref:`sec:CMake`.


.. _sec:GNUMake:

Building with GNU Make
======================

Here, we walk you through compiling an incflo executable.
Decide whether you want to build the executable in 2-D vs 3-D,
with or without EB (embedded boundaries / cut cells), and choose the corresponding
directory: ``test_2d``, ``test_3d``, ``test_no_eb_2d``, or ``test_no_eb_3d``.
For this example, let's say you want to build in 2-D with cut cells.
Then to build the code:

#. ``cd`` to the desired build directory; e.g.

   .. code:: shell

             cd incflo/test_2d


#. Edit the ``GNUmakefile``:

   Set ``AMREX_HOME`` to be the path to the directory where you have put amrex, and ``AMREX_HYDRO_HOME``
   to be the path to the directory where you have put AMReX-Hydro. These default to ``../../amrex`` and
   ``../../AMReX-Hydro``, i.e. to amrex, AMReX-Hydro and incflo having been cloned into the same directory.
   NOTE: when setting ``AMREX_HOME`` in the ``GNUmakefile``, be aware that ``~`` does not expand, so ``AMREX_HOME=~/amrex/`` will yield an error.

   Alternatively, the paths can be set up as environment variables, ``AMREX_HOME`` and ``AMREX_HYDRO_HOME``, on your machine to point to the path names where you have put AMReX and AMReX-Hydro. For example, if you are using the bash shell, you can add this to your ``.bashrc`` as:

   .. highlight:: bash

   ::

      export AMREX_HOME=/path/to/amrex
      export AMREX_HYDRO_HOME=/path/to/AMReX-Hydro

   alternatively, in tcsh one can set

   .. highlight:: tcsh

   ::

      setenv AMREX_HOME /path/to/amrex
      setenv AMREX_HYDRO_HOME /path/to/AMReX-Hydro


   Common options that you can set in the ``GNUmakefile`` include (the default values shown are
   those set in the test directories' ``GNUmakefile``, or AMReX's defaults for options not set there)

   +-----------------+------------------------------+------------------+-------------+
   | Option name     | Description                  | Possible values  | Default     |
   |                 |                              |                  | value       |
   +=================+==============================+==================+=============+
   | COMP            | Compiler (gnu or intel)      | gnu / intel      | gnu         |
   +-----------------+------------------------------+------------------+-------------+
   | USE_MPI         | Whether to enable MPI        | TRUE / FALSE     | TRUE        |
   +-----------------+------------------------------+------------------+-------------+
   | USE_OMP         | Whether to enable OpenMP     | TRUE / FALSE     | FALSE       |
   +-----------------+------------------------------+------------------+-------------+
   | USE_CUDA        | Whether to enable CUDA       | TRUE / FALSE     | FALSE       |
   +-----------------+------------------------------+------------------+-------------+
   | DEBUG           | Whether to use DEBUG mode    | TRUE / FALSE     | FALSE       |
   +-----------------+------------------------------+------------------+-------------+
   | PROFILE         | Include profiling info       | TRUE / FALSE     | FALSE       |
   +-----------------+------------------------------+------------------+-------------+
   | TINY_PROFILE    | Include tiny profiling info  | TRUE / FALSE     | FALSE       |
   +-----------------+------------------------------+------------------+-------------+
   | COMM_PROFILE    | Include comm profiling info  | TRUE / FALSE     | FALSE       |
   +-----------------+------------------------------+------------------+-------------+
   | TRACE_PROFILE   | Include trace profiling info | TRUE / FALSE     | FALSE       |
   +-----------------+------------------------------+------------------+-------------+

   .. note::
      **Do not set both USE_OMP and USE_CUDA to true.**

   Additional options and information on using other compilers can be found in the
   AMReX documentation at
   https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html .


#. Make the executable:

   Now type

   .. code:: shell

      make

   The name of the resulting executable (generated by the make system) encodes several of the build characteristics, including dimensionality of the problem, compiler name, and whether MPI and/or OpenMP were linked with the executable.
   Thus, several different build configurations may coexist simultaneously in a problem folder.
   For example, the default build in ``incflo/test_2d`` will look
   like ``incflo2d.gnu.MPI.EB.ex``, indicating that this is a 2-d version of the code, made with
   ``COMP=gnu``, ``USE_MPI=TRUE`` and support for embedded boundaries (EB).


.. _sec:CMake:

Building with CMake
===================

CMake build is a two-steps process. First ``cmake`` is invoked to create
configuration files and makefiles in a chosen directory (``builddir``).
Next, the actual build is performed by invoking ``make`` from within ``builddir``.
The CMake build process is summarized as follows:

.. highlight:: console

::

    mkdir /path/to/builddir
    cd    /path/to/builddir
    cmake [options] -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] /path/to/incflo
    make

In the above snippet, ``[options]`` indicates one or more options for the
customization of the build, as described later in this section.
If the option ``CMAKE_BUILD_TYPE`` is omitted,
``CMAKE_BUILD_TYPE=Release`` is assumed.

There are two modes to build incflo with CMake:

* **SUPERBUILD (recommended):** AMReX and AMReX-Hydro are built as part
  of the incflo build process. This method is strongly encouraged as it
  ensures that the configuration options are consistent.

* **STANDALONE:** incflo source code is built separately and linked to an existing
  AMReX installation (AMReX-Hydro is still built as part of incflo).
  This is ideal for continuous integration servers (CI)
  and regression testing applications. AMReX library version and configuration options
  must meet incflo requirements.



.. note::
   **incflo requires CMake 3.20 or higher.**

.. _sec:build:superbuild:

SUPERBUILD Instructions (recommended)
-------------------------------------

In this mode, incflo CMake inherits AMReX CMake targets and configuration options, as well as AMReX-Hydro CMake targets; that is, incflo, AMReX and AMReX-Hydro are configured and built as a single entity.

First, make sure you've downloaded all the source codes (see :ref:`sec:Download`).

Next, build the incflo executable:

   .. code-block:: console

             > cd incflo
             > mkdir build
             > cd build
             > cmake -DAMREX_HOME=/path/to/amrex/source/directory [amrex options] -DAMREX_HYDRO_HOME=/path/to/AMReX-Hydro/source/directory [incflo options] -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] ..
             > make -j


where ``AMREX_HOME`` is a CMake variable pointing to the top-level source directory of the AMReX distribution you downloaded earlier; ``[amrex options]`` is a list of any of the AMReX configuration options listed in the `AMReX user guide - Building with CMake <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-with-cmake>`_; ``AMREX_HYDRO_HOME`` is a CMake variable pointing to the top-level source directory of the AMReX-Hydro distribution you downloaded earlier; ``[incflo options]`` is a list of the incflo
configuration options listed in the table below; and ``..`` is the relative path pointing to the top-level source directory of incflo or this example.
incflo sets AMReX's dimension, MPI, OpenMP, GPU backend, EB and particle options from the
corresponding incflo options (e.g. ``-DINCFLO_DIM=2``, not ``-DAMReX_SPACEDIM=2``), and always builds
AMReX in double precision with its linear solvers.


.. _sec:build:standalone:

**STANDALONE** instructions
---------------------------------------------------------------------

For a stand alone build, incflo CMake will look for an existing
AMReX installation on the system and link the incflo binaries against it.

Building AMReX
~~~~~~~~~~~~~~~~~~~

First, make sure you've downloaded all the source codes (see :ref:`sec:Download`).

Next, configure, build and install AMReX as follows:

.. code:: shell

          > cd /path/to/amrex
          > mkdir build
          > cd build
          > cmake -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel]  [amrex options] -DCMAKE_INSTALL_PREFIX:PATH=/absolute/path/to/installdir ..
          > make install

``[amrex options]`` is a list of any of the AMReX configuration options listed in the
`AMReX user guide - Building with CMake <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-with-cmake>`_.
AMReX must be built with the configuration that incflo will use: the same dimension
(``AMReX_SPACEDIM``, 3 by default), MPI, OpenMP and GPU backend, EB (``-DAMReX_EB=YES``) and particles
if incflo uses them, double precision and the linear solvers (AMReX's defaults). CMake checks this
when configuring incflo.


Building incflo
~~~~~~~~~~~~~~~~~


Build incflo:

.. code:: shell

    > cd /path/to/incflo
    > mkdir build
    > cd build
    > cmake -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] [incflo options] -DAMReX_ROOT=/absolute/path/to/amrex/installdir -DAMREX_HYDRO_HOME=/path/to/AMReX-Hydro/source/directory ..
    > make -j


Passing ``-DAMReX_ROOT=/absolute/path/to/amrex/installdir`` instructs CMake to search
``/absolute/path/to/amrex/installdir`` before searching system paths
for an available AMReX installation.
``AMReX_ROOT`` can also be set as an environmental variable instead of passing it as a command line option.
``AMREX_HYDRO_HOME`` is a CMake variable pointing to the top-level source directory of the AMReX-Hydro
distribution you downloaded earlier; it is required in both modes, since AMReX-Hydro is always built
as part of incflo.


``[incflo options]`` indicates any of the configuration option listed in the table below.

+-----------------+------------------------------+------------------+-------------+
| Option name     | Description                  | Possible values  | Default     |
|                 |                              |                  | value       |
+=================+==============================+==================+=============+
| CMAKE\_CXX\     | User-defined C++ flags       | valid C++        | None        |
| _FLAGS          |                              | compiler flags   |             |
+-----------------+------------------------------+------------------+-------------+
| CMAKE\_CUDA\    | User-defined CUDA flags      | valid CUDA       | None        |
| _FLAGS          |                              | compiler flags   |             |
+-----------------+------------------------------+------------------+-------------+
| INCFLO\_DIM     | Dimensionality of the build  | 2/3              | 3           |
|                 |                              |                  |             |
+-----------------+------------------------------+------------------+-------------+
| INCFLO\_MPI     | Enable build with MPI        | no/yes           | yes         |
|                 |                              |                  |             |
+-----------------+------------------------------+------------------+-------------+
| INCFLO\_OMP     | Enable build with OpenMP     | no/yes           | no          |
|                 |                              |                  |             |
+-----------------+------------------------------+------------------+-------------+
| INCFLO\_CUDA    | Enable build with CUDA       | no/yes           | no          |
|                 |                              |                  |             |
+-----------------+------------------------------+------------------+-------------+
| INCFLO\_HIP     | Enable build with HIP        | no/yes           | no          |
|                 |                              |                  |             |
+-----------------+------------------------------+------------------+-------------+
| INCFLO\_SYCL    | Enable build with SYCL       | no/yes           | no          |
|                 |                              |                  |             |
+-----------------+------------------------------+------------------+-------------+
| INCFLO\_EB      | Build Embedded Boundary      | no/yes           | no          |
|                 | support                      |                  |             |
+-----------------+------------------------------+------------------+-------------+
| INCFLO\_CSG     | Build CSG geometry support   | no/yes           | no          |
|                 | (requires INCFLO_EB)         |                  |             |
+-----------------+------------------------------+------------------+-------------+
| INCFLO\_HYPRE   | Enable HYPRE support         | no/yes           | no          |
|                 |                              |                  |             |
+-----------------+------------------------------+------------------+-------------+
| INCFLO\_FPE     | Build with Floating-Point    | no/yes           | no          |
|                 | Exceptions checks            |                  |             |
+-----------------+------------------------------+------------------+-------------+
| INCFLO\_\       | Build Tracer Particles       | no/yes           | no          |
| PARTICLES       | support                      |                  |             |
+-----------------+------------------------------+------------------+-------------+




Few more notes on building incflo
-----------------------------------

The system default C++ compiler can be overwritten as follows:

.. code:: shell

    > cmake -DCMAKE_CXX_COMPILER=<c++-compiler>  [options]  ..

When building on a platform that uses the ``module`` utility, use either
the above command (with full path to the compilers) or the following:

.. code:: shell

    > cmake -DCMAKE_CXX_COMPILER=CC [options] ..

incflo uses the same compiler flags used to build AMReX, unless
``CMAKE_CXX_FLAGS`` is explicitly provided, or
the environmental variables ``CXXFLAGS`` is set.



Running the code
================
incflo takes an input file as its first command-line argument.  The file may
contain a set of parameter definitions that will override defaults set in the code.
Many example inputs files can be found in the ``incflo/test_XXX`` directories.
For example, to run the Taylor-Green vortices example in ``incflo/test_2d`` on 4 MPI processes with the
executable built there with GNU Make (see :ref:`sec:GNUMake`), type:

.. code:: shell

          > cd incflo/test_2d
          > mpiexec -n 4 ./incflo2d.gnu.MPI.EB.ex benchmark.taylor_green_vortices

With CMake, the executable is ``incflo.ex`` in the build directory.


For more information on inputs options, see the :ref:`Chap:Inputs` section.

incflo can generate subfolders in the current folder
named ``plt00000``, ``plt00010``, etc, and ``chk00000``,
``chk00010``, etc. These are called plotfiles and checkpoint
files. The plotfiles are used for visualization of derived fields; the checkpoint
files are used for restarting the code. For more information see the AMReX Documentation on
:ref:`amrex:Chap:IO`.


Visualizing the Results
=======================

Options for visualizing the data are discussed in the AMReX Documentation on :ref:`amrex:Chap:Visualization`.
