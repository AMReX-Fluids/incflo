Fluid Variables
===============

   +-----------------------+---------------------------------------------------------------------+
   | Variable              | Definition                                                          |
   +=======================+=====================================================================+
   | :math:`\rho`          | Fluid density                                                       |
   +-----------------------+---------------------------------------------------------------------+
   | :math:`U`             | Fluid velocity                                                      |
   +-----------------------+---------------------------------------------------------------------+
   | :math:`p`             | Pressure                                                            |
   +-----------------------+---------------------------------------------------------------------+
   | :math:`\eta`          | Viscosity (``incflo.mu`` for a Newtonian fluid)                     |
   +-----------------------+---------------------------------------------------------------------+
   | :math:`\tau`          | Viscous stress tensor, :math:`\eta (\nabla U + (\nabla U)^T)`       |
   +-----------------------+---------------------------------------------------------------------+
   | :math:`\mu_s`         | scalar diffusivity                                                  |
   +-----------------------+---------------------------------------------------------------------+
   | :math:`\mu_T`         | thermal conductivity                                                |
   +-----------------------+---------------------------------------------------------------------+
   | :math:`c_p`           | Specific heat (``incflo.cp``)                                       |
   +-----------------------+---------------------------------------------------------------------+
   | :math:`{\bf g}`       | Gravitational acceleration                                          |
   +-----------------------+---------------------------------------------------------------------+
   | :math:`{\bf H}_U`     | :math:`= (H_x , H_y , H_z )`, External Forces                       |
   +-----------------------+---------------------------------------------------------------------+
   | :math:`H_s`           | External sources                                                    |
   +-----------------------+---------------------------------------------------------------------+
   | :math:`H_T`           | External heat sources                                               |
   +-----------------------+---------------------------------------------------------------------+

Fluid Equations
===============

Conservation of fluid mass:

.. math:: \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho U)  = 0

Conservation of fluid momentum:

.. math:: \frac{ \partial (\rho U)}{\partial t}
   + \nabla \cdot (\rho U U) + \nabla p = \nabla \cdot \tau  + {\bf H}_U

By default, incflo solves the equivalent convective form,
:math:`\rho ( \partial U / \partial t + U \cdot \nabla U ) + \nabla p = \nabla \cdot \tau + {\bf H}_U`;
``incflo.advect_momentum = true`` selects the conservative form above.
With ``incflo.use_tensor_solve = false``, the velocity components are diffused separately, with
:math:`\nabla \cdot (\eta \nabla U)` in place of :math:`\nabla \cdot \tau` (the two are equal for a
constant viscosity).

Velocity constraint:

.. math:: \nabla \cdot U = 0

or, with ``incflo.use_divu_source = 1``, the constraint

.. math:: \nabla \cdot U = S

for a given source :math:`S`, as in low Mach number models. Mass is still conserved, so
the density changes along particle paths as :math:`D\rho/Dt = -\rho S`. With
``incflo.constant_density = 1`` the density is not evolved, so mass is not conserved:
the constraint adds mass at the rate :math:`\rho_0 S` per unit volume. That is not a low
Mach number model, and incflo warns about it, but it has exact solutions that are useful
for testing. :math:`S` comes from ``incflo::compute_divu_source``, which by default
uses an analytic :math:`S` defined for the chosen ``incflo.probtype``. It is virtual,
so a class derived from ``incflo`` can compute :math:`S` from the state instead.
In a closed domain (only walls and periodic boundaries) the
projections can only enforce :math:`\nabla \cdot U = S - \bar{S}`, where :math:`\bar{S}` is
the mean of :math:`S`; incflo warns when :math:`\bar{S} \ne 0`. That is a uniform change
of the thermodynamic pressure, which incflo does not evolve.

Tracer(s): for conservative,

.. math:: \frac{\partial \rho s}{\partial t} + \nabla \cdot (\rho U s)  =  \nabla \cdot \mu_s \nabla s + \rho H_s

or, for non-conservative (see ``incflo.trac_is_conservative``),

.. math:: \frac{\partial s}{\partial t} + U \cdot \nabla s  =  \nabla \cdot \mu_s \nabla s + H_s

Optional temperature:

.. math:: \rho c_p \left( \frac{\partial T}{\partial t} + U \cdot \nabla T \right)  =  \nabla \cdot \mu_T \nabla T + H_T

By default, :math:`H_s = 0`,  :math:`H_T = 0`, and :math:`{\bf H}_U = {\bf 0}`.
If gravity is set (``incflo.gravity``), then :math:`{\bf H}_U = \rho {\bf g}`, except for the
Boussinesq problems (``incflo.probtype`` = 11, 111, 112, 113 and 500), which use
:math:`{\bf H}_U = \rho s_0 {\bf g}`, where :math:`s_0` is the first tracer.

.. _sec:BackgroundPressure:

Background pressure
-------------------

Except for the Boussinesq problems, incflo writes the pressure as :math:`p = p_0 + p'`, where the
background pressure :math:`p_0` has a constant gradient, and solves for :math:`p'`:

- in a direction with gravity, :math:`\nabla p_0 = \rho_0 {\bf g}`, the hydrostatic pressure gradient
  for the density :math:`\rho_0` (``incflo.ro_0``);
- in the direction in which ``incflo.delp`` is non-zero, :math:`\nabla p_0 = - \mathrm{delp} / L`,
  where :math:`L` is the length of the domain in that direction. This pressure drop drives the flow,
  e.g. in a periodic channel;
- in a direction with a pressure inflow at one end and a pressure outflow at the other,
  :math:`\nabla p_0 = (p_{hi} - p_{lo}) / L`, where :math:`p_{lo}` and :math:`p_{hi}` are the
  boundary pressures.

Only one direction can have a ``delp`` or a pressure inflow/outflow pair, and gravity cannot be in
that direction. The pressure written to plotfiles (``p``) is :math:`p'`, while the pressure gradient
(``gpx``, ``gpy`` and ``gpz``) includes :math:`\nabla p_0`.
