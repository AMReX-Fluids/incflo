.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _Chap:EB:

Embedded Boundaries
===================

incflo follows AMReX's approach to embedded boundaries (EB), which is described in the
`AMReX EB documentation <https://amrex-codes.github.io/amrex/docs_html/EB_Chapter.html>`_. 
By default, incflo uses AMReX's constructive solid geometry framework defined in the namespace :cpp:`amrex::EB2`.
Alternatively, the constructive solid geometry can also be created using OpenSCAD's CSG file format by installing
the :cpp:`csg-eb` library. To use this option, incflo must be built with the flag
:cpp:`USE_CSG=TRUE` and :cpp:`CSGEB_HOME` must be set to where the library was installed.
See `MFIX's CSG-EB repository <https://mfix.netl.doe.gov/gitlab/exa/csg-eb>`_ for more details about this format.

incflo provides several options of embedded boundary geometries. The inputs parameter ``incflo.geometry = XXX``
determines which geometry is selected by :cpp:`incflo::MakeEBGeometry()` within :cpp:`incflo/src/embedded_boundaries`.
The procedure to create your own EB geometry is described in the AMReX documentation on :ref:`amrex:sec:EB:ebinit`.
As discussed in the AMReX documentation, note that when constructing the EB, we must specify a
maxium coarsening level (``max_crse_level``):

.. highlight:: c++

::

   EB2::Build(gshop, geom[lev], required_crse_lev, max_crse_level);


This specifies to which level of coarseness the EB is still defined. It might not be
immediately obvious, but the multigrid solver (used in the fluid solve) also
depends indirectly on this parameters. Choosing a value of ``max_crse_level`` that is too small might restrict
how many levels the MLMG solver can use, and therefore give slightly different answers in the fluid solve.
