.. _Chap:InputsPlotfiles:

Plotfiles and Other Output
==========================

The following inputs must be preceded by "amr" and control frequency and naming of plotfile generation as well
as whether the EB geometry should be written out.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| plot_int            | Frequency of plotfile output;                                         |    Int      | -1        |
|                     | if -1 then no plotfiles will be written at this frequency             |             |           |
| plot_per_exact      | Time period of plotfile output (exact); will modify dt                |    Real     | -1        |
|                     | if -1 then no plotfiles will be written at this frequency             |             |           |
| plot_per_approx     | Time period of plotfile output (approximate); does not modify dt      |    Real     | -1        |
|                     | if -1 then no plotfiles will be written at this frequency             |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plotfile_on_restart | Should we write a plotfile when we restart (only used if plot_int>0)  |   Bool      | False     |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plot_file           | Prefix to use for plotfile output                                     |  String     | plt       |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| write_eb_surface    | Should we write out the EB geometry in vtp format                     |   Bool      | False     |
|                     | If true, it will only be written once,after initialization or restart |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+

The following inputs must be preceded by "amr" and control what variables will be written in plotfiles.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| plt_regtest         | Save all variables to plot file (overrides all other IO flags)        |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_vel             | Save fluid velocity data to plot file                                 |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_p               | Save fluid pressure to plot file                                      |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_ro              | Save fluid density to plot file                                       |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_mu              | Save fluid viscosity to plot file                                     |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_gradp           | Save gradient of pressure filed to plot file                          |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_vort            | Save vorticity to plot file                                           |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
