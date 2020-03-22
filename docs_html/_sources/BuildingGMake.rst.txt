Building MFiX-Exa with gmake
============================

Building MFiX-Exa with gmake 

If you want to use gmake to build MFiX_Exa, you will need to have already
cloned amrex into a local directory:

.. code:: shell

    > git clone https://github.com/amrex-codes/amrex

If you want to run MFIX-Exa using the SIMPLE algorithm or the projection
method with normal velocity components defined on faces, then the easiest
way to build it is:

.. code:: shell

    > git clone http://mfix.netl.doe.gov/gitlab/exa/mfix.git
    > cd mfix/exec

(If you want to run MFIX-Exa using the projection method with cell-centered
velocity components, replace exec above with exec_cc)

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

Then type

.. code:: shell

    > make
