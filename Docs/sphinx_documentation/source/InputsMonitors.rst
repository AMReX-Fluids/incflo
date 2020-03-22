.. _Chap:InputsMonitors:

Monitors
========

The following inputs must be preceded by "amr" and control whether to compute
monitors, i.e., spatial averages, and how often to output the results.
n is the number of monitors implicitly defined by the size of avg_region_x_w.

+------------------+-----------------------------------------------------------------------+-------------+-----------+
|                  | Description                                                           |   Type      | Default   |
+==================+=======================================================================+=============+===========+
| avg_int          | Interval, in number of CFD dt's, to write output                      |  Int        | -1        |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| avg_file         | Base file name which is appended with the data type (vel_p, p_g,      |  String     | avg_region|
|                  | ep_g or vel_g), the number of this type of averaging,  and the .csv   |             |           |
|                  | file extension                                                        |             |           |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| avg_vel_p        | Average and save particle velocity (if 1)                             |  n*Int      | 0         |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| avg_p_g          | Average and save fluid-phase pressure (if 1)                          |  n*Int      | 0         |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| avg_ep_g         | Average and save fluid-phase volume fraction (if 1)                   |  n*Int      | 0         |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| avg_vel_g        | Average and save fluid-phase velocity (if 1)                          |  n*Int      | 0         |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| avg_region_x_w   | Lower bound of averaging region in x-direction                        |  n*Real     | None      |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| avg_region_x_e   | Upper bound of averaging region in x-direction                        |  n*Real     | None      |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| avg_region_y_s   | Lower bound of averaging region in y-direction                        |  n*Real     | None      |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| avg_region_y_n   | Upper bound of averaging region in y-direction                        |  n*Real     | None      |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| avg_region_z_b   | Lower bound of averaging region in z-direction                        |  n*Real     | None      |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| avg_region_z_t   | Upper bound of averaging region in z-direction                        |  n*Real     | None      |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
