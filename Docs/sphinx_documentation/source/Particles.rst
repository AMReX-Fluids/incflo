.. role:: cpp(code)
   :language: c++

.. _Chap:Particles:

Tracer Particles
================

incflo provides Lagrangian particle support built on AMReX's :cpp:`ParticleContainer` machinery, which is
described in the `AMReX particle documentation <https://amrex-codes.github.io/amrex/docs_html/Particle_Chapter.html>`_.
The particle framework is used to seed the flow with massless tracer particles that advect with the fluid
velocity. These particles can be used to visualize flow structures, to map a mass density back onto the
mesh for post-processing, and as a refinement criterion for adaptive mesh refinement.

Particles are organized into species, each of which is an :cpp:`incflo_PC` object held by the
:cpp:`ParticleData` container in :cpp:`incflo/src/particles`. At present incflo defines a single
species, ``tracer_particles``, whose run-time parameters are read from the ``tracer_particles.*`` input
prefix. The particle module is compiled in with the :cpp:`INCFLO_USE_PARTICLES` preprocessor flag.

Enabling tracer particles
--------------------------

Tracer particles are disabled by default. They are activated by setting the following parameter, which
must be preceded by "incflo":

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| use_tracer_particles | Enables the use of tracer particles that advect with the flow.        |    Int      | 0            |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

Tracer particle parameters
---------------------------

The following inputs must be preceded by "tracer_particles" and control how the particles are initialized
and advected:

+------------------------------+-----------------------------------------------------------------------+-----------------+--------------+
|                              | Description                                                           |      Type       | Default      |
+==============================+=======================================================================+=================+==============+
| initial_distribution_type    | Type of particle initialization; currently the only supported         |     String      | box          |
|                              | option is "box", i.e., a uniform distribution of particles in a box.  |                 |              |
+------------------------------+-----------------------------------------------------------------------+-----------------+--------------+
| particle_box_lo              | Low corner of the box in which particles are initially placed;        |  Vector<Real>   | ProbLo       |
|                              | defaults to the low corner of the computational domain.               |                 |              |
+------------------------------+-----------------------------------------------------------------------+-----------------+--------------+
| particle_box_hi              | High corner of the box in which particles are initially placed;       |  Vector<Real>   | ProbHi       |
|                              | defaults to the high corner of the computational domain.              |                 |              |
+------------------------------+-----------------------------------------------------------------------+-----------------+--------------+
| place_randomly_in_cells      | If true, particles are placed at random positions within each cell;   |      Bool       | true         |
|                              | if false, a fixed offset is used, which is useful for regression      |                 |              |
|                              | testing.                                                              |                 |              |
+------------------------------+-----------------------------------------------------------------------+-----------------+--------------+
| initial_particles_per_cell   | Initial number of particles to place in each grid cell.               |      Int        | 1            |
+------------------------------+-----------------------------------------------------------------------+-----------------+--------------+
| advect_with_flow             | If true, particles are advected with the flow velocity using a        |      Bool       | true         |
|                              | second-order midpoint method; the default is true for tracer          |                 |              |
|                              | particles.                                                            |                 |              |
+------------------------------+-----------------------------------------------------------------------+-----------------+--------------+

A complete example of the tracer particle inputs is

::

   incflo.use_tracer_particles = 1

   tracer_particles.initial_distribution_type = box
   tracer_particles.initial_particles_per_cell = 4
   tracer_particles.particle_box_lo = 0.0 0.0 0.0
   tracer_particles.particle_box_hi = 1.0 1.0 1.0
   tracer_particles.advect_with_flow = true

Cylinder constraint
-------------------

The particle framework includes an optional cylinder constraint that removes particles that leave a
specified cylindrical region. It is a demonstration of how to keep particles inside a region of interest
(e.g., flow in a cylinder) and is controlled by the following inputs, which must be preceded by
"cylinder":

+---------------+-----------------------------------------------------------------------+-----------------+--------------+
|               | Description                                                           |      Type       | Default      |
+===============+=======================================================================+=================+==============+
| radius        | Radius of the cylinder inside which particles are kept; if set to a   |      Real       | none         |
|               | value greater than zero, particles outside the cylinder are removed   |                 |              |
|               | at initialization and after each advection step.                      |                 |              |
+---------------+-----------------------------------------------------------------------+-----------------+--------------+
| center        | Coordinates of a point on the cylinder axis.                          |  Vector<Real>   | none         |
+---------------+-----------------------------------------------------------------------+-----------------+--------------+
| direction     | Axis along which the cylinder is oriented: 0 = x, 1 = y, 2 = z.       |      Int        | none         |
+---------------+-----------------------------------------------------------------------+-----------------+--------------+

Initialization
--------------

When tracer particles are enabled, they are initialized as soon as a level is created from scratch in
:cpp:`incflo::MakeNewLevelFromScratch()`. The only initialization supported at present is a uniform
distribution in a box (``initial_distribution_type = box``), which is implemented in
:cpp:`incflo_PC::initializeParticlesUniformDistributionInBox()` in :cpp:`incflo/src/particles/incflo_PCInit.cpp`:

- Each cell of the box is assigned ``initial_particles_per_cell`` particles placed at random positions
  within the cell, unless ``place_randomly_in_cells = false``, in which case a fixed offset is used.
- If incflo is built with embedded boundary support, particles are only placed in cells with a non-zero
  volume fraction.
- All particles are initialized with zero velocity and a mass of :math:`10^{-6}`.
- If the cylinder inputs are present, particles lying outside the cylinder are removed immediately after
  initialization.

Note that, as currently implemented, particles are placed only in cells whose indices satisfy
:math:`i \bmod 2 = 0` and :math:`j \bmod 2 = 1`, so only one quarter of the cells in the box receive
particles.

Advection
---------

At every time step, the tracer particles are advanced by :cpp:`incflo::evolveTracerParticles()`, which is
called from the predictor step (unless the advection scheme is :cpp:`MOL`) and from the corrector step of
the projection algorithm. Each level is advanced with the same time step :math:`\Delta t` as the fluid.

For each particle species with ``advect_with_flow = true``, :cpp:`incflo_PC::AdvectWithFlow()` advances
the particles with a second-order midpoint method using the MAC face velocities:

#. In a first pass, the flow velocity at the particle position is interpolated from the MAC face
   velocities and each particle is moved by half a time step to a midpoint position; the original
   position is saved.
#. In a second pass, the velocity is interpolated at the midpoint position and the particle is moved the
   full time step from its saved original position.

The velocity interpolation is performed with :cpp:`mac_interpolate` in :cpp:`AMReX_TracerParticle_mod_K.H`.
On non-periodic boundaries the particles are reflected back into the domain and the sign of the
corresponding velocity component is reversed. After advection the particle velocity attribute is set to
the interpolated flow velocity, and :cpp:`Redistribute()` is called to move particles between MPI ranks
and levels as needed. If the cylinder constraint is enabled, particles that leave the cylinder are
removed by setting their ids to -1.

Output
------

Particle data are written to checkpoint files and to plotfiles using AMReX particle I/O. In addition, the
following mesh variables can be derived from the particles:

- A particle count per cell (computed with :cpp:`incflo_PC::Increment()`), which is included in plotfiles
  as ``particle_count`` and is used in :cpp:`incflo::ErrorEst()` as a refinement criterion: when
  ``incflo.refine_particles`` is true (the default), cells containing at least one particle are tagged
  for refinement.
- A mass density computed by depositing the particle mass onto the mesh with linear interpolation in
  :cpp:`incflo_PC::massDensity()`, which appears in plotfiles as the variable
  ``tracer_particles_mass_density``.
