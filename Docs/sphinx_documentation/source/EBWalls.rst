.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran


.. _sec:EB-basics:


Constructing Embedded Boundaries in MFiX-Exa
============================================

MFiX uses AMReX's constructive solid geometry framework defined in the namespace
:cpp:`amrex::EB2`. See the `AMReX EB documentation`_ for more details. These are
defined in ``src/eb/mfix_eb.cpp``. A the function :cpp:`mfix::make_eb_geometry`
(also defined in ``src/eb/mfix_eb.cpp``) selects :cpp:one of the following
geometries depending on the value of the :cpp:``mfix.geometry`` setting in the
``inputs`` file.

+------------------------------+----------------------+-------------------------+
|   Description                |  ``mfix.geometry``   |   Implementation Satus  |
+==============================+======================+=========================+
| Planar walls from mfix.dat   | Don't specify        | Fully implemented       |
| (on by default)              |                      |                         |
+------------------------------+----------------------+-------------------------+
| Box (up to six walls)        | ``box``              | Fully implemented       |
+------------------------------+----------------------+-------------------------+
| Cylinder                     | ``cylinder``         | Fully implemented       |
+------------------------------+----------------------+-------------------------+
| Hopper                       | ``hopper``           | Fully implemented       |
+------------------------------+----------------------+-------------------------+
| Cyclone                      | ``cyclone``          | Fully implemented       |
+------------------------------+----------------------+-------------------------+
| General                      | ``general``          | Partially implemented   |
|                              | cf. note 1 below     |                         |
+------------------------------+----------------------+-------------------------+
| Hourglass                    | ``hourglass``        | Deprecated              |
|                              | cf. note 1           | cf. note 2 below        |
+------------------------------+----------------------+-------------------------+
| CLR (chemical looping        | ``clr``              | Deprecated              |
| reactor)                     | cf. note 1           | cf. note 2              |
+------------------------------+----------------------+-------------------------+
| CLR Riser                    | ``clr_riser``        | Deprecated              |
|                              | cf. note 1           | cf. note 2              |
+------------------------------+----------------------+-------------------------+

1. Older (legacy) alternative settings are:

   +-----------------------------+-------------------------------+
   | Value of  ``mfix.geometry`` |  Alternative                  |
   +=============================+===============================+
   | ``general``                 | ``mfix.use_poly2 = true``     |
   |                             | or ``mfix.use_walls = true``  |
   +-----------------------------+-------------------------------+
   | ``hourglass``               | ``mfix.hourglass = true``     |
   +-----------------------------+-------------------------------+
   | ``clr``                     | ``mfix.clr = true``           |
   +-----------------------------+-------------------------------+
   | ``clr_riser``               | ``mfix.clr_riser = true``     |
   +-----------------------------+-------------------------------+

2. These geometries where not ported from AMReX's old :cpp:``EB`` system to the
   new :cpp:``EB2``.

Also note that planar boundary conditions can be specified in the ``mfix.dat``
file. Even if the user does not specify an ``mfix.geometry`` in the ``inputs``,
any no-slip or free-slip boundary conditions are expressed as EB walls.

There are also two parameters (specified in the ``inputs`` file) that influence
the level-set creation:

+-------------------------------+------------------------------------------------+
|  Parameter                    |  Description                                   |
+===============================+================================================+
| ``amr.max_level``             | If greater than 0, MFiX operates in multi-level|
|                               | mode. The level-set grids follow all other     |
|                               | grids. If equal to 1, the level-set has two    |
|                               | levels (one additional level with higher       |
|                               | refinement).                                   |
+-------------------------------+------------------------------------------------+
| ``mfix.levelset__refinement`` | If ``amr.max_level > 1`` this parameter is     |
|                               | ignored. Otherwise it sets the maximum         |
|                               | refinement of the level-set                    |
+-------------------------------+------------------------------------------------+


How MFiX-Exa Constructs the EB Geometry
---------------------------------------

Once a geometry is selected by :cpp:`mfix::make_eb_geometry`, the procedure is
the same for (almost) all geometries. Also see the `AMReX geometry
documentation`_ for information on how to construct new geometries:

1. Construct an implicit function representing the geometry (using the language
   of constructive solid geometry). For example

.. highlight:: c++

::

   EB2::CylinderIF my_cyl(radius, height, direction, center, inside);
   auto gshop_cyl = EB2::makeShop(my_cyl);

2. (Optional) Construct the implicit function representing the EB seen by the
   particles. This might deviate from the "standard" EB depending on the
   specific application. In all standard cases used by mfix, this step is
   omitted.

3. Call :cpp:`mfix::build_eb_levels(gshop)` this function builds the EB levels
   and fills the implicit function :cpp:`MultiFab` (the later being used to
   construct the level-set function). Note that this also makes the particle EB
   levels point to the fluid eb levels.

4. (Optional, do this if you did 2.) Call
   :cpp:`mfix::build_particle_eb_levels(gshop_part)` **last**. This will update
   the particle EB levels.


MFiX's EB Data Structures
-------------------------

The :cpp:`mfix` class stores the following EB data:

.. highlight:: c++

::

   //! EB levels representing fluid boundary conditions
   Vector<const EB2::Level *> eb_levels;
   //! EB levels representing particle boundary conditions (same as
   //! `mfix::eb_levels` but might include additional walls).
   Vector<const EB2::Level *> particle_eb_levels;

   //! EB factory that lives on the fluid grids
   Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > ebfactory;
   //! EB factory that lives on the particle grids
   Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > particle_ebfactory;

As discussed in the previous sub-section, the difference between
:cpp:`mfix::eb_levels` and :cpp:`mfix::particle_eb_levels` enables the user to
specify a modfied EB geometry for particles only. Whereas the fluid sees the EB
geometry in :cpp:`mfix::eb_levels`. If no addition particle EB geometry is
specified (point 4 in the previous section), then
:cpp:`mfix::particle_eb_levels` points to :cpp:`mfix::eb_levels`.

In the same spirit, the :cpp:`mfix::ebfactory` is constructed over the fluid
grid and using the fluid EB levels, whereas :cpp:`mfix::particle_ebfactory` is
constructed over the particle grid using the particle EB levels.


A note about constructing EB Levels
-----------------------------------

MFiX-Exa builds EB levels in :cpp:`mfix::build_eb_levels` (via
:cpp:`LSCore<F>::BuildEBLevel`)

.. highlight:: c++

::

   EB2::Build(gshop, geom[lev], required_crse_lev, max_crse_level);
   const EB2::IndexSpace & ebis = EB2::IndexSpace::top();


When building an EB level, the maximum coarsening level (:cpp:`int
max_crse_level`) and the required coarsening level (:cpp:`int
required_crse_lev`) need to be specified. The reason for this is that we need to
specify to which level of coarseness the EB is still defined. It might not be
immediately obvious, but the Poisson solver (used in the fluid solve) also
depends indirectly on these parameters. Thus changing these during EB level
creation might restrict how many levels the MLMG solver can use, and therefore
give slightly different answers in the fluid solve.



Local Mesh Refinement at Walls
==============================

MFiX-Exa has the capability of locally refining the computational grid near EBs.
This is done by tagging (in :cpp:`mfix::ErrorEst`) any cells with volume
fraction between 0 and 1. To enable local mesh refinement, set ``amr.max_level``
to a value greater than 1. Note that the parameter ``mfix.levelset__refinement``
is ignored on all cases except when ``amr.max_level = 1``.


MFiX-Exa Initialization Process
-------------------------------

Since MFiX requires the volume fraction when building grids (because this is
needed by :cpp:`mfix::ErrorEst`), the EB geometries need to be built before
calling :cpp:`mfix::Init`. The recommended procedure therefore is

.. highlight:: c++

::

   // Default constructor (geom[lev] is defined here)
   mfix my_mfix;

   // Initialize internals from ParamParse database
   my_mfix.InitParams(solve_fluid, solve_dem, call_udf);

   // Initialize memory for data-array internals
   my_mfix.ResizeArrays();

   // Construct EB (must be done _before_ mfix::Init)
   my_mfix.make_eb_geometry();

   // Initialize derived internals. Grids are create here.
   my_mfix.Init(dt, time);

   // Create EB factories on new grids
   my_mfix.make_eb_factories();

   if (solve_dem)
   {
       // Fill level-sets on each level (must be done _after_ mfix::Init)
       my_mfix.fill_eb_levelsets();
   }

   // Finish constructing levels
   my_mfix.InitLevelData(dt,time);

   // Regrid (ensure all MultiFabs are on their correct grids)
   my_mfix.Regrid();


Also note that mfix defines boundary conditions in Fortran also (via the
mfix.dat). Since these are potentially needed to build EB walls,
:cpp:`mfix::make_eb_geometry` also calls :cpp:`mfix_set_bc_type`.
 
The grids for each level are build in the :cpp:`mfix::Init` by invoking the
initialization functions inherited from :cpp:`amrex::AmrCore`.

.. highlight:: c++

::

   // This tells the AmrMesh class not to iterate when creating the initial
   // grid hierarchy
   SetIterateToFalse();

   // This tells the Cluster routine to use the new chopping routine which
   // rejects cuts if they don't improve the efficiency
   SetUseNewChop();

   // This Builds the new Grids
   InitFromScratch(0.);



The Level-Set Function
======================

MFiX-Exa uses a level-set function to resolve particle-wall collisions. See the
`AMReX Level-Set documentation`_ for more details. The level-set function is
stored on the nodal :cpp:`Vector<std::unique_ptr<MultiFab>> mfix::level_sets`.
The level-set data is always stored on the particle grids. Depending on the
input ``amr.max_level`` The level-set can be in one of two modes:

1. MFiX-Exa is running in single-level mode (:cpp:`nlev == 1`). Then
   :cpp:`mfix::level_sets[0]` will be at the same resolution as the fluid
   (except that it is stored on the particle grid). Even though :cpp:`nlev == 1`,
   there is a second level, :cpp:`level_sets[1]`. This level is the same as
   :cpp:`level_sets[0]` but refined by :cpp:`mfix::levelset__refinement`. This
   way the level-set always has the appropriate resolution to resolve structures
   in the EB, even if the fluid is defined on a fairly coarse grid.

2. MFiX-Exa is running in multi-level mode (:cpp:`nlev > 1`). The the parameter
   :cpp:`mfix::levelset__refinement` is ignored. :cpp:`mfix::level_sets` then
   follows the rest of MFiX, i.e. it is defined on the particle grids on all
   levels.

The level-set is used in two places:

1. The function :cpp:`MFIXParticleContainer::EvolveParticles` interpolates the
   level-set onto each particle's position in order to resolve collisions with
   the EBs. If :cpp:`nlev == 1`, :cpp:`level_sets[1]` is used to evolve the
   particle positions. Otherwise :cpp:`level_sets[lev]` is used for each level.

2. The fluid-particle coupling can sometimes rely on neighbor stencils where one
   or more cell is covered by an EB. In order to avoid values that do not
   conform with the boundary conditions, the fluid velocity is reconstructed in
   those cells. The algorithm relies on the level-set, and uses
   :cpp:`level_sets[lev]` on each level.


Special Cases Involving Level-Sets
----------------------------------

The level-set function is filled by the `mfix::fill_eb_levelsets()` function.
There are two special cases involving level-sets:

1. Mass-Inflow boundary conditions are not given EB walls. However, we don't
   want particles to fall out of a MI either, so at the very end of the
   `mfix::fill_eb_levelsets()` function we call `mfix::intersect_ls_walls()`.
   This performs an intersection operation with the level-set representing a
   wall at each MI.

2. Box geometries and regular geometries are comprised entirely out of planar
   surfaces. Therefore the levelset is not construction out of an EB factory (as
   would be the case for all other geometries). But out of an intersection with
   all planar surfaces. This has the advantage of correctly describing corners.

.. _AMReX EB documentation: https://amrex-codes.github.io/amrex/docs_html/EB_Chapter.html
.. _AMReX Level-Set documentation: https://amrex-codes.github.io/amrex/docs_html/EB.html#level-sets
.. _AMReX geometry documentation: https://amrex-codes.github.io/amrex/docs_html/EB.html#initializing-the-geometric-database 
