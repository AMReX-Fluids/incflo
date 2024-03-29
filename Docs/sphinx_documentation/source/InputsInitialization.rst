.. _Chap:InputsInitialization:

Initialization
==============

The following inputs must be preceded by "amr" and determine how we initialize a calculation:

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| restart              | If set, then restart from this file rather than from scratch          |  String     |   None       |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+


The following inputs must be preceded by "incflo" and determine how we initialize a calculation:

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| do_initial_proj      | Should we do the initial projection?                                  |    Bool     |  True        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| initial_iterations   | How many pressure iterations before starting the first timestep       |  Int        |    3         |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
