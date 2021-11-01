Building incflo with gmake
============================

Building incflo with gmake

If you want to use gmake to build incflo, you will need to have already
cloned amrex into a local directory:

.. code:: shell

    > git clone https://github.com/amrex-codes/amrex

Then,

.. code:: shell

    > git clone http://github.com/amrex-codes/incflo.git
    > cd incflo

Now decide whether you want to build the executable in 2-D vs 3-D,
with EB (embedded boundaries / cut cells) or not.  For this example,
let's say you want to build in 2-D with cut cells.  Then

.. code:: shell

    > git clone http://github.com/amrex-codes/incflo.git
    > cd test_2d

(The other directories are test_3d, test_no_eb, and test_no_eb_2d.)

Edit the GNUmakefile to set AMREX_HOME to be the path to the directory
where you have put amrex.  Other options that you can set include

+-----------------+------------------------------+------------------+-------------+
| Option name     | Description                  | Possible values  | Default     |
|                 |                              |                  | value       |
+=================+==============================+==================+=============+
| COMP            | Compiler (gnu or intel)      | gnu / intel      | None        |
+-----------------+------------------------------+------------------+-------------+
| USE_MPI         | Whether to enable MPI        | TRUE / FALSE     | FALSE       |
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

Then type

.. code:: shell

    > make
