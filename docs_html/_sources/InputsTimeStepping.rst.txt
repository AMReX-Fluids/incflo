.. role:: cpp(code)
   :language: c++

.. _sec:InputsTimeStepping:

Time Stepping
=============

The first three inputs below do not take a prefix.  Note that if the first two are both specified, both criteria
are used and the simulation will stop when the first criterion is hit, i.e. when either the number of steps
reaches max_step or time reaches stop_time.

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| max_step             | Maximum number of time steps to take                                  |    Int      |  -1          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| stop_time            | Maximum time to reach                                                 |    Real     | -1.0         |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| steady_state         | Is this a steady state problem? Not yet implemented: setting it to    |    Bool     |   False      |
|                      | true aborts when the inputs are read                                  |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

The inputs below must be preceded by "incflo."

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| fixed_dt             | Value of fixed dt if > 0                                              |    Real     |   -1.        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| cfl                  | CFL constraint (dt < cfl * dx / u) if fixed_dt not > 0;               |    Real     |   0.5        |
|                      | at most 0.5 for MOL and 1.0 for Godunov and BDS                       |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| dt_change_max        | Maximum factor by which dt can increase from one step to the next;    |    Real     |   1.1        |
|                      | must be between 1.0 and 1.1                                           |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| init_shrink          | reduce the initial time step by this safety factor                    |    Real     |   0.1        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

Setting the Time Step
---------------------

  * If you want to fix the dt, simply set :cpp:`incflo.fixed_dt = XXX` and the fluid time
    step will always be that number.

  * If you want to let the code determine the appropriate time step using the advective CFL
    condition, then set :cpp:`incflo.cfl = 0.7` for example, and (in the absence of forcing
    terms and of explicit diffusion, which further restrict it) the fluid time step will
    be computed to be dt = 0.7 * dx / max(vel).

  * Note that the cfl defaults to 0.5 so it does not have to be set in the inputs file. If neither
    :cpp:`incflo.cfl` nor :cpp:`incflo.fixed_dt` is set, then default value of cfl will be used.
    If :cpp:`incflo.fixed_dt` is set, then it will override the cfl option whether
    :cpp:`incflo.cfl` is set or not.
