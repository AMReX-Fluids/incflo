Building incflo with gmake
============================

Building incflo with gmake 

If you want to use gmake to build incflo, you will need to have already
cloned amrex into a local directory:

.. code:: shell

    > git clone https://github.com/amrex-codes/amrex

Then

.. code:: shell

    > git clone http://github.com/amrex-codes/incflo.git
    > cd incflo/test

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
