Particles on GPUs
==========================

The particle components of MFIX-Exa are a natural candidate for offloading to the GPU. 
The particle kernels are compute-intensive and can in principle be processed asynchronously with parts of the fluid advance.
The core components of the particle method in MFIX-Exa are:

- Neighbor List Construction
- Particle-Particle Collisions
- Particle-Wall Collisions

Of these operations, the neighbor list construction requires the most care. 
A neighbor list is a pre-computed list of all the neighbors a given particle can interact with over the next *n* timesteps. 
Neighbor lists are usually constructed by binning the particles by an interaction distance, 
and then performing the N\ :sup:`2` distance check only on the particles in neighboring bins. In detail, the CPU version of the neighbor list algorithm is as follows:

- For each tile on each level, loop over the particles, identifying the bin it belongs to.
- Add the particle to a linked-list for the cell that `owns` it.
- For each cell, loop over all the particles, and then loop over all potential collisions partners in the neighboring cells.
- If a collision partner is close enough, add it to that particle's neighbor list.

To port this algorithm to the GPU, we use the parallel algorithms library Thrust, distributed as part of the CUDA Toolkit. Thrust provides parallel sorting, searching, and prefix summing algorithms that are particularly useful in porting particle algorithms. To construct the neighbor list on the GPU, we follow the basic approach used by Canaba, a product of the Particle Co-Design Center:

- Sort the particles on each grid by bin, using a parallel counting sort. We use Thrust's `exclusive\_scan` function to implement the prefix sum phase of the sort, and hand-coded kernels for the rest. This step does not actually involving rearranging the particle data - rather, we compute a permutation that would put the particles in order without actually reordering them.
- Once the particles are sorted by bin, we can loop over the particles in neighboring bins. We make two passes over the particles. First, we launch a kernel to count the number of collision partners for each particle.
- Then, we sum these numbers and allocate space for our neighbor list.
- Finally, we make a another pass over the particles, putting them into to list at the appropriate place.

Note that we build a \emph{full} neighbor list, meaning that if particle $i$ appears in particle $j$'s list, then particle $j$ also appears in particle $i$'s list. This simplifies the force-computation step when using these lists, since the forces and torques for a given particle can be updated without atomics.

The final on-grid neighbor list data structure consists of two arrays. First, we have the neighbor list itself, stored as a big, 1D array of particle indices. Then, we have an `offsets` array that stores, for each particle, where in the neighbor list array to look. The details of this data structure have been hidden inside an iterator, so that user code can look like:

.. code-block:: c

       // now we loop over the neighbor list and compute the forces                                                                                                                                                
        AMREX_FOR_1D ( np, i,
        {
            ParticleType& p1 = pstruct[i];
            p1.rdata(PIdx::ax) = 0.0;
            p1.rdata(PIdx::ay) = 0.0;
            p1.rdata(PIdx::az) = 0.0;

            for (const auto& p2 : nbor_data.getNeighbors(i))
            {
                Real dx = p1.pos(0) - p2.pos(0);
                Real dy = p1.pos(1) - p2.pos(1);
                Real dz = p1.pos(2) - p2.pos(2);

                ...
            }

Note that, because of our use of managed memory to store the particle data and the neighbor list, the above code will work when compiled for either CPU or GPU.

The above algorithm deals with constructing a neighbor list for the particles on a single grid. When domain decomposition is used, one must also make copies of particles on adjacent grids, potentially performing the necessary MPI communication for grids associated with other processes. The routines `fillNeighbors`, which computes which particles needed to be ghosted to which grid, and `updateNeighbors`, which copies up-to-date data for particles that have already been ghosted, have also been offloaded to the GPU, using techniques similar to AMReX's `Redistribute` routine. The important thing for users is that calling these functions does not trigger copying data off the GPU.

Once the neighbor list has been constructed, collisions with both particles and walls can easily be processed. 

We have created a GPU branch of MFIX that is capable of running with GPU support. As of this writing, the following operations in MFIX have been offloaded:

- Neighbor particles / neighbor list construction
- Particle-particle collisions
- Particle-wall collisions
- PIC Deposition (used in putting the drag force and solids volume fraction on the grid)



