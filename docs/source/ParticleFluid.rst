.. role:: cpp(code)
   :language: c++


.. _sec:particle-fluid:

Particle/Fluid Interactions
===========================

A standard :cpp:`mfix_level::Evolve` step looks like this (for brevity, we are
omitting some details here):

.. highlight:: c++

::

    // Calculate particle volume fraction:
    mfix_calc_volume_fraction(lev, sum_vol)

    // Eolve fluid
    if ( use_proj_method )
        EvolveFluidProjection(lev, nstep, dt, ...)
    else
        EvolveFluidSimple(lev, nstep, dt, ...)

    // Apply fluid drag force on particles
    mfix_calc_drag_particle(lev)

    // Move particles (using forces acting on particles)
    pc -> EvolveParticles(lev, nstep, dt ...)

Here, ``lev`` represents the AMR level, and ``pc`` is a pointer to a
:cpp:`MFIXParticleContainer` instance.
