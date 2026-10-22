.. _Chap:InputsPlotfiles:

Plotfiles and Other Output
==========================

The following inputs must be preceded by "amr" and control frequency and naming of plotfile generation.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| plot_int            | Frequency of plotfile output;                                         |    Int      | -1        |
|                     | if -1 then no plotfiles will be written at this frequency             |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plot_per_exact      | Time period of plotfile output (exact); will modify dt                |    Real     | -1        |
|                     | if -1 then no plotfiles will be written at this frequency             |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plot_per_approx     | Time period of plotfile output (approximate); does not modify dt      |    Real     | -1        |
|                     | if -1 then no plotfiles will be written at this frequency.            |             |           |
|                     | Only one of plot_int, plot_per_exact and plot_per_approx may be > 0   |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plotfile_on_restart | Should we write a plotfile when we restart                            |   Bool      | False     |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plot_file           | Prefix to use for plotfile output                                     |  String     | plt       |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+

The following input must be preceded by "incflo" and controls whether the EB geometry is written out.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| write_eb_surface    | Should we write out the EB geometry in vtp format (3D EB builds only) |   Bool      | False     |
|                     | If true, it will only be written once, after initialization or        |             |           |
|                     | restart.                                                              |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+

The following inputs must be preceded by "amr" and control what variables will be written in plotfiles.
By default, incflo plotfiles contain the velocity, pressure gradient, density, tracer, velocity magnitude, vorticity,
and if using EB, volume fraction. With ``incflo.use_divu_source = 1`` they also contain ``divu_source``, the source
:math:`S` in the divergence constraint :math:`\nabla \cdot U = S`.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| plotVariables       | Space separated list of variable names. **Takes precedence over all** |  String     | see above |
|                     | **other plotfile options**. Plotfile will contain only this list,     |             |           |
|                     | e.g. "amr.plotVariables = velx density" will result in plotfiles      |             |           |
|                     | containing only the x-component of the velocity and the density.      |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_ccse_regtest    | Collection of variables for regression test plotfiles (will get       |    Int      | 0         |
|                     | overwritten by plotVariables or plt_* flags)                          |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_*               | The plt_* flags below are deprecated in favor of plotVariables.       |    Int      |           |
|                     | Each one adds (1) or removes (0) one variable from the default list.  |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_velx            | Save x-velocity to plot file                                          |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_vely            | Save y-velocity to plot file                                          |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_velz            | Save z-velocity to plot file                                          |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_p               | Save pressure to plot file                                            |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_rho             | Save density to plot file                                             |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_tracer          | Save tracer to plot file                                              |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_eta             | Save viscosity to plot file                                           |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_gpx             | Save dp/dx to plot file                                               |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_gpy             | Save dp/dy to plot file                                               |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_gpz             | Save dp/dz to plot file                                               |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_magvel          | Save velocity magnitude to plot file                                  |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_vort            | Save vorticity to plot file                                           |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_forcing         | Save forcing term for velocity to plot file                           |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_strainrate      | Save strain rate to plot file                                         |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_divu            | Not implemented: requesting divu is refused when the inputs are read  |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_vfrac           | Save EB volume fraction to plot file (EB builds only)                 |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+

If ``incflo.use_temperature = 1``, the temperature is also added to the default list of plotfile variables.

incflo provides the option to write a user defined selection of variables to a "smallplotfile" at an interval that is
different than that of the full plotfile.
The following inputs must be preceded by "amr" and control frequency and naming of smallplotfile generation.

+--------------------------+-----------------------------------------------------------------------+-------------+-----------+
|                          | Description                                                           |   Type      | Default   |
+==========================+=======================================================================+=============+===========+
| smallplot_int            | Frequency of smallplotfile output;                                    |    Int      | -1        |
|                          | if -1 then no smallplotfiles will be written at this frequency        |             |           |
+--------------------------+-----------------------------------------------------------------------+-------------+-----------+
| smallplot_per_approx     | Time period of smallplotfile output (approximate); does not modify dt |    Real     | -1        |
|                          | if -1 then no smallplotfiles will be written at this frequency.       |             |           |
|                          | Only one of smallplot_int and smallplot_per_approx may be > 0         |             |           |
+--------------------------+-----------------------------------------------------------------------+-------------+-----------+
| smallplotfile_on_restart | Write smallplotfile when restarting                                   |    Bool     | False     |
+--------------------------+-----------------------------------------------------------------------+-------------+-----------+
| smallplot_file           | Prefix to use for smallplotfile output                                |  String     | smallplt  |
+--------------------------+-----------------------------------------------------------------------+-------------+-----------+
| smallplotVariables       | Space separated list of variable names.                               |  String     | none      |
+--------------------------+-----------------------------------------------------------------------+-------------+-----------+
