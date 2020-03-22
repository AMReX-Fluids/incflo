.. _Chap:Particles:

Particles
=========

MFiX-Exa supports solid particles which are immersed in the fluid. In the
future particles can be implemented using the immersed-boundary method (and
hence these IBM-particles can support non-spherical shapes, and special
hydrodynamical and chemical interactions). But for now, particles are
spherical, with approximate fluid-particle interactions. The main consequence
of this approximation is that particles should not be bigger than the fluid
grid resolution.

Particles interact with container walls (implemented using AMReX's
Embedded-Boundary functionality) using a level-set constructed from the signed
closest distance function at the beginning of the simulation.

Refer to the subsequent sections for details on the implementation of the
particle-particle, particle-fluid, and particle-wall interactions.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   ParticleBasics
   ParticleFluid
   ParticleWalls
   ParticlesOnGpus
