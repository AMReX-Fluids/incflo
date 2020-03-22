.. role:: cpp(code)
   :language: c++


.. _sec:particle-basics:

Particle Basics
===============

In MFiX-Exa, particles are managed by the `MFIXParticleContainer
<https://amrex-codes.github.io/MFIX-Exa/doxygen/class_m_f_i_x_particle_container.html>`_
class.  This class is derived from AMReX's :cpp:`NeighborParticleContainer`
and handles all of the particle data. 
:cpp:`MFIXParticleContainer` provides the functions for solving the
particle dynamics (based on particle-particle, particle-fluid, and
particle-wall forces)


Particle Dynamics
-----------------

During the DES steps, particle positions are updated using the
`MFIXParticleContainer::EvolveParticles
<https://amrex-codes.github.io/MFIX-Exa/doxygen/class_m_f_i_x_particle_container.html#a158f3f5fa11c262ad6a9909b40a5cd13>`_
method. It's structure is:

.. code-block:: c++
   :linenos:
   :emphasize-lines: 29-39, 42-45

    // Set time-step size (subdt) and number (numbsebsteps) for the DES steps
    des_init_time_loop( &time, &dt, &nsubsteps, &subdt, &subdt_io );

    // Temporary storage of forces and torques
    std::map<PairIndex, Vector<Real>> tow;
    std::map<PairIndex, Vector<Real>> fc;
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        tow[index] = Vector<Real>();
        fc[index] = Vector<Real>();
    }

    while (n < nsubsteps) // step over number of substeps (DES part)
    {
        // Neighborlist house-keeping
        if (n % 25 == 0) {
            clearNeighbors(lev);
            Redistribute();
            fillNeighbors(lev);
            buildNeighborList(lev,sort_neighbor_list);
        } else {
            updateNeighbors(lev);
        }

        // Itterate over particles
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            // Forces and torques due to particle-wall collisions
            calc_wall_collisions_ls(particles, & ntot, & nrp,
                                    tow[index].dataPtr(), fc[index].dataPtr(), & subdt,
                                    ...
                                    );

            calc_particle_collisions(particles                     , &nrp,
                                     neighbors[index].dataPtr()    , &size_ng,
                                     neighbor_list[index].dataPtr(), &size_nl,
                                     tow[index].dataPtr(), fc[index].dataPtr(),
                                     &subdt, &ncoll);
        }

        // Move particles based on velocities and forces
        des_time_loop(&nrp     , particles,
                      &ntot, tow[index].dataPtr(), fc[index].dataPtr(), &subdt,
                      &xlen, &ylen, &zlen, &stime, &n);
    }

where the highlighted lines are responsible for moving the particles.
