Building incflo with CMake
============================

CMake build is a two-step process. First ``cmake`` is invoked to create
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

There are two modes to build incflo with cmake:

o **SUPERBUILD (recommended):** AMReX is downloaded and built as part
of the incflo build process. This method is strongly encouraged as it
ensures that the configuration options for incflo and AMReX are consistent.

o **STANDALONE:** incflo source code is built separately and linked to an existing
AMReX installation. This is ideal for continuous integration severs (CI)
and regression testing applications. AMReX library version and configuration options
must meet incflo requirements.



.. note::
   **incflo requires CMake 3.14 or higher.**

.. _sec:build:superbuild:

SUPERBUILD Instructions (recommended)
-------------------------------------

By default incflo CMake looks for an existing AMReX installation on the system. If none is found, it falls back to **SUPERBUILD** mode.
In this mode, incflo CMake inherents AMReX CMake targets and configuration options, that is, incflo and AMReX are configured and built as a single entity

Assuming no valid AMReX installation is present on the target system, and ``AMReX_ROOT`` is not set (see :ref:`sec:build:standalone`), the following code will build incflo in **SUPERBUILD** mode:

.. code:: shell

    > git clone http://github.com/amrex-codees/incflo.git
    > cd incfo
    > mkdir build
    > cd build
    > cmake [amrex options] -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] ..
    > make -j

``[amrex options]`` in the snippet above is a list of any of the AMReX configuration options listed in the `AMReX user guide <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-with-cmake>`_


Working with the AMReX submodule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SUPERBUILD** mode relies on a git submodule to checkout the AMReX git repository.
If the AMReX submodule is not initialized, **SUPERBUILD** mode will initialize it and checkout 
the AMReX commit the submodule is pointing to.
Instead, if the AMReX submodule has already been manually initialized and a custom commit has been checked out,
**SUPERBUILD** mode will use that commit. For incflo development or testing, you may need to build with a different
branch or release of AMReX.

The ``subprojects/amrex`` directory is a Git repo, so use all standard Git
commands. Either ``cd subprojects/amrex`` to run Git commands in the ``amrex``
directory, or use ``git -C subprojects/amrex`` in the incflo repo. For
instance, to build with the ``my-amrex-branch`` branch of the AMReX repo:

.. code:: shell

    > git -C subprojects/amrex checkout my-amrex-branch
    > git status
    ...
    modified:   subprojects/amrex (new commits)

The branch ``my-amrex-branch`` will then be used when building incflo.

To revert to the default version of the AMReX submodule, run ``git submodule
update``:

.. code:: shell

    > git submodule update
    > git status
    ...
    nothing to commit, working tree clean

You can edit, commit, pull, and push AMReX changes from ``subprojects/amrex``.
AMReX development is outside the scope of this document. Run ``git status`` in
the top-level incflo repo to see whether the AMReX submodule has new commits,
modified files, or untracked files.

The CMake variables ``AMREX_GIT_COMMIT_DEVELOP`` and
``AMREX_GIT_COMMIT_DEVELOP`` have been removed. Instead, to update the AMReX
submodule referenced by incflo:

.. code:: shell

    > git -C subprojects/amrex checkout UPDATED_AMREX_COMMIT_SHA1
    > git add subprojects/amrex
    > git commit -m 'Updating AMReX version'

This will only update the AMReX SHA-1 referenced by incflo. Uncommitted AMReX
changes and untracked AMReX files under ``subprojects/amrex`` are not added by
``git add subprojects/amrex``. (To commit to the AMReX repo, change directories
to ``subprojects/amrex`` and run Git commands there, before ``git add
subprojects/amrex``.)

.. note::

    Only update the AMReX submodule reference in coordination with the other
    incflo developers!

   
.. _sec:build:standalone:

**STANDALONE** instructions
---------------------------------------------------------------------

Building AMReX 
~~~~~~~~~~~~~~~~~~~

Clone AMReX from the official Git repository and checkout the
*development* branch:

.. code:: shell

    > git clone https://github.com/AMReX-Codes/amrex.git
    > cd amrex
    > git checkout development

Next, configure, build and install AMReX as follows:

.. code:: shell

    > mkdir build
    > cd build
    > cmake -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] -DENABLE_AMRDATA=yes -DENABLE_EB=yes [other amrex options] -DCMAKE_INSTALL_PREFIX:PATH=/absolute/path/to/installdir ..
    > make install

The options  **ENABLE\_AMRDATA=yes**, and **ENABLE\_EB=yes** are required by incflo. ``[other amrex options]`` in the snippet above refers to any other AMReX configuration option in addition to the required ones. Please refer to the `AMReX user guide <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-with-cmake>`_ for more details on building AMReX with CMake.


Building incflo
~~~~~~~~~~~~~~~~~

Clone and build incflo:

.. code:: shell

    > git clone https://github.com/amrex-codes/incflo.git
    > mkdir build
    > cd build
    > cmake -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] [incflo options] -DAMReX_ROOT=/absolute/path/to/amrex/installdir ..
    > make -j


Passing ``-DAMReX_ROOT=/absolute/path/to/amrex/installdir`` instructs CMake to search
``/absolute/path/to/amrex/installdir`` before searching system paths
for an available AMReX installation.
``AMReX_ROOT`` can also be set as an environmental variable instead of passing it as a command line option.


``[incflo options]`` indicates any of the configuration option listed in the table below.

+-----------------+------------------------------+------------------+-------------+
| Option name     | Description                  | Possible values  | Default     |
|                 |                              |                  | value       |
+=================+==============================+==================+=============+
| CMAKE\_Fortran\ | User-defined Fortran flags   | valid Fortran    | None        |
| _FLAGS          |                              | compiler flags   |             |
+-----------------+------------------------------+------------------+-------------+
| CMAKE\_CXX\     | User-defined C++ flags       | valid C++        | None        |
| _FLAGS          |                              | compiler flags   |             |
+-----------------+------------------------------+------------------+-------------+
| CMAKE\_CUDA\    | User-defined CUDA flags      | valid CUDA       | None        |
| _FLAGS          |                              | compiler flags   |             |
+-----------------+------------------------------+------------------+-------------+
| ENABLE\_MPI     | Enable build with MPI        | no/yes           | yes         |
|                 |                              |                  |             |
+-----------------+------------------------------+------------------+-------------+
| ENABLE\_OMP     | Enable build with OpenMP     | no/yes           | no          |
|                 |                              |                  |             |
+-----------------+------------------------------+------------------+-------------+
| ENABLE\_CUDA    | Enable build with CUDA       | no/yes           | no          |
|                 |                              |                  |             |
+-----------------+------------------------------+------------------+-------------+
| ENABLE\_HYPRE   | Enable HYPRE support         | no/yes           | no          |
|                 |                              |                  |             |
+-----------------+------------------------------+------------------+-------------+
| ENABLE\_FPE     | Build with Floating-Point    | no/yes           | no          |
|                 | Exceptions checks            |                  |             |
+-----------------+------------------------------+------------------+-------------+




Few more notes on building incflo
-----------------------------------

The system defaults compilers can be overwritten as follows:

.. code:: shell

    > cmake -DCMAKE_CXX_COMPILER=<c++-compiler> -DCMAKE_Fortran_COMPILER=<f90-compiler> [options]  ..

When building on a platform that uses the ``module`` utility, use either
the above command (with full path to the compilers) or the following:

.. code:: shell

    > cmake -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn [options] ..

incflo uses the same compiler flags used to build AMReX, unless
``CMAKE_Fortran_FLAGS``/``CMAKE_CXX_FLAGS`` is explicitly provided, or
the environmental variables ``FFLAGS``/``CXXFLAGS`` are set.


Building incflo for Cori (NERSC)
-----------------------------------

Standard build 
~~~~~~~~~~~~~~~~~~~

For the Cori cluster at NERSC, you first need to load/unload modules required to build incflo.

.. code:: shell

    > module unload altd
    > module unload darshan
    > module load cmake/3.14.0

The default options for Cori are the **Haswell** architecture and **Intel** compiler, if you want to compile with the **Knight's Landing (KNL)** architecture:

.. code:: shell

    > module swap craype-haswell craype-mic-knl

Or use the **GNU** compiler:

.. code:: shell

    > module swap PrgEnv-intel PrgEnv-gnu

Now incflo can be built following the :ref:`sec:build:superbuild`.

.. note::

    The load/unload modules options could be saved in the `~/.bash_profile.ext` 


GPU build 
~~~~~~~~~~~~~~~~~~~

To compile on the GPU nodes in Cori, you first need to purge your modules, most of which won't work on the GPU nodes

.. code:: shell

    > module purge

Then, you need to load the following modules:

.. code:: shell

    > module load modules esslurm gcc cuda openmpi/3.1.0-ucx cmake/3.14.0

Currently, you need to use OpenMPI; mvapich2 seems not to work.

Then, you need to use slurm to request access to a GPU node:

.. code:: shell

    > salloc -N 1 -t 02:00:00 -c 80 -C gpu -A m1759 --gres=gpu:8 --exclusive

This reservers an entire GPU node for your job. Note that you can’t cross-compile for the GPU nodes - you have to log on to one and then build your software.

Finally, navigate to the base of the incflo repository and compile in GPU mode:

.. code:: shell

    > cd incflo
    > mdkir build
    > cd build
    > cmake -DENABLE_CUDA=yes -DCUDA_ARCH=Volta -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran ..
    > make -j

For more information about GPU nodes in Cori -- `<https://docs-dev.nersc.gov/cgpu/>`_

Building incflo for Summit (OLCF)
-----------------------------------

For the Summit cluster at OLCF, you first need to load/unload modules required to build incflo.

.. code:: shell

    > module unload xalt
    > module unload darshan
    > module load gcc
    > module load cmake/3.14.0

Now incflo can be built following the :ref:`sec:build:superbuild`.

To build incflo for GPUs, you need to load cuda module:

.. code:: shell

    > module load cuda/10.1.105

To compile for GPUs:

.. code:: shell

    > cd incflo
    > mdkir build
    > cd build
    > cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DENABLE_CUDA=yes 
    > make -j

An example of a *submission_script* for GPUs can be found in the repo ``incflo/tests/GPU_test/script.sh``.
For more information about Summit cluster: `<https://www.olcf.ornl.gov/for-users/system-user-guides/summit/>`_ 
