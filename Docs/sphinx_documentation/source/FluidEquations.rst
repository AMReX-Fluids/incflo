Fluid Variables
===============

   +-----------------------+--------------------------------------------------+
   | Variable              | Definition                                       |
   +=======================+==================================================+
   | :math:`\rho`          | Fluid density                                    |
   +-----------------------+--------------------------------------------------+
   | :math:`U`             | Fluid velocity                                   |
   +-----------------------+--------------------------------------------------+
   | :math:`\tau`          | Viscous stress tensor                            |
   +-----------------------+--------------------------------------------------+
   | :math:`\mu_s`         | scalar diffusivity                               |
   +-----------------------+--------------------------------------------------+
   | :math:`{\bf g}`       | Gravitational acceleration                       |
   +-----------------------+--------------------------------------------------+
   | :math:`{\bf H}_U`     | :math:`= (H_x , H_y , H_z )`, External Forces    |
   +-----------------------+--------------------------------------------------+
   | :math:`H_s`           | External sources                                 |
   +-----------------------+--------------------------------------------------+

Fluid Equations
===============

Conservation of fluid mass:

.. math:: \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho U)  = 0

Conservation of fluid momentum:

.. math:: \frac{ \partial (\rho U)}{\partial t}
   + \nabla \cdot (\rho U U) + \nabla p = \nabla \cdot \tau  + {\bf H}_U

Incompressibility constraint:

.. math:: \nabla \cdot U = 0

Tracer(s): for conservative,

.. math:: \frac{\partial \rho s}{\partial t} + \nabla \cdot (\rho U s)  =  \nabla \cdot \mu_s \nabla s + \rho H_s

or, for non-conservative,

.. math:: \frac{\partial s}{\partial t} + U \cdot \nabla s  =  \nabla \cdot \mu_s \nabla s + H_s

By default, :math:`H_s = 0` and :math:`{\bf H}_U = {\bf 0}`.
If gravity is set during runtime, then :math:`{\bf H}_U` defaults to :math:`\rho {\bf g}`.
