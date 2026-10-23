Time Step
=========

incflo has the option to treat advection using either the Method of Lines (MOL) or a Godunov
approach.
Here we discuss the basic workflow involved in a time step for each of these approaches.
For details on how the convective terms are constructed, see the AMReX-Hydro documentation on :ref:`hydro:schemes`.
For details on the projections, see the AMReX-Hydro documentation on :ref:`hydro:projections`.

The updates below are written for the conservative form of the momentum equation, which is used
with ``incflo.advect_momentum = true``. By default, incflo updates the velocity instead: the
advective term is :math:`U^{MAC} \cdot \nabla U`, and the pressure gradient, viscous and forcing terms
are divided by the density. :math:`{\bf H}_U^{n+1/2}` is the forcing term evaluated with
:math:`\rho^{n+1/2}`, and the density in the projections is :math:`\rho^{n+1/2}`.

MOL
~~~~~~~~~~~~~~~~

MOL requires predictor-corrector methodology to achieve second order accuracy.
In the predictor

-  Define :math:`U^{MAC,n}`, the face-centered (staggered) MAC velocity which is used for advection, using :math:`U^n`

-  Define an approximation to the new-time state, :math:`(\rho U)^{\ast}` by setting

.. math:: (\rho U)^{\ast} =& (\rho U)^n -
           \Delta t \left( \nabla \cdot (\rho U^{MAC} U)^n + \nabla {p}^{n-1/2} \right) \\ &+
           \Delta t \left( \nabla \cdot \tau^\ast + {\bf H}_U^{n+1/2} \right)

for implicit diffusion, which is the default. For Crank-Nicolson diffusion
(``incflo.diffusion_type = 1``), :math:`\Delta t \; \nabla \cdot \tau^\ast` is replaced by
:math:`\frac{\Delta t}{2} ( \nabla \cdot \tau^n + \nabla \cdot \tau^\ast )`, and for explicit
diffusion (``incflo.diffusion_type = 0``) by :math:`\Delta t \; \nabla \cdot \tau^n`.

-  Project :math:`U^{\ast}` by solving

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi = \nabla \cdot \left( \frac{1}{\Delta t}
          U^{\ast}+ \frac{1}{\rho} \nabla {p}^{n-1/2} \right)

then defining

.. math:: U^{\ast \ast} = U^{\ast} + \frac{\Delta t}{\rho} \left( \nabla {p}^{n-1/2} - \nabla \phi \right)

and

.. math:: {p}^{n+1/2, \ast} = \phi


In the corrector

-  Define :math:`U^{MAC,\ast \ast}` at the "new" time using :math:`U^{\ast \ast}`

-  Define a new approximation to the new-time state, :math:`(\rho U)^{\ast \ast \ast}` by setting

.. math:: (\rho U)^{\ast \ast \ast} =& (\rho U)^n - \frac{\Delta t}{2} \left( \nabla \cdot (\rho U^{MAC} U)^n + \nabla \cdot (\rho U^{MAC} U)^{\ast \ast}\right) \\
          &+ \Delta t \; \nabla \cdot \tau^{\ast \ast \ast} + \Delta t \left( - \nabla {p}^{n+1/2,\ast} + {\bf H}_U^{n+1/2} \right)

for implicit diffusion. For Crank-Nicolson diffusion, :math:`\Delta t \; \nabla \cdot \tau^{\ast \ast \ast}`
is replaced by :math:`\frac{\Delta t}{2} ( \nabla \cdot \tau^n + \nabla \cdot \tau^{\ast \ast \ast} )`, and
for explicit diffusion by :math:`\frac{\Delta t}{2} ( \nabla \cdot \tau^n + \nabla \cdot \tau^{\ast \ast} )`.

-  Project :math:`U^{\ast \ast \ast}` by solving

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi = \nabla \cdot \left( \frac{1}{\Delta t} U^{\ast \ast \ast} + \frac{1}{\rho} \nabla {p}^{n+1/2,\ast} \right)

then defining

.. math:: U^{n+1} = U^{\ast \ast \ast} + \frac{\Delta t}{\rho} \left( \nabla {p}^{n+1/2,\ast} - \nabla \phi \right)

and

.. math:: {p}^{n+1/2} = \phi

Godunov Methods
~~~~~~~~~~~~~~~~~~~~

When we use a time-centered Godunov approach (i.e. the ``Godunov`` or ``BDS`` options),
we no longer need the predictor and corrector steps.

-  Define the time-centered face-centered (staggered) MAC velocity which is used for advection: :math:`U^{MAC,n+1/2}`

-  Define the new-time density, :math:`\rho^{n+1} = \rho^n - \Delta t \nabla \cdot (\rho^{n+1/2,pred} U^{MAC,n+1/2})`

-  Define an approximation to the new-time state, :math:`(\rho U)^{\ast}` by setting

   .. math:: (\rho^{n+1} U^{\ast}) &= (\rho^n U^n) -
           \Delta t \left( \nabla \cdot (\rho U^{MAC} U)^{n+1/2} + \nabla {p}^{n-1/2} \right)  \\
           &+ \Delta t \; \nabla \cdot \tau^\ast
            + \Delta t \; {\bf H}_U^{n+1/2}

   for implicit diffusion, which is the current default. For Crank-Nicolson diffusion
   (``incflo.diffusion_type = 1``), :math:`\Delta t \; \nabla \cdot \tau^\ast` is replaced by
   :math:`\frac{\Delta t}{2} ( \nabla \cdot \tau^n + \nabla \cdot \tau^\ast )`, and for explicit
   diffusion (``incflo.diffusion_type = 0``) by :math:`\Delta t \; \nabla \cdot \tau^n`.

-  Project :math:`U^{\ast}` by solving

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi = \nabla \cdot \left( \frac{1}{\Delta t}
          U^{\ast}+ \frac{1}{\rho} \nabla {p}^{n-1/2} \right)

then defining

.. math:: U^{n+1} = U^{\ast} + \frac{\Delta t}{\rho} \left( \nabla {p}^{n-1/2} - \nabla \phi \right)

and

.. math:: {p}^{n+1/2} = \phi

With ``incflo.use_divu_source = 1``, the MAC projection enforces
:math:`\nabla \cdot U^{MAC,n+1/2} = S^{n+1/2}`, and :math:`S^{n+1} / \Delta t` is subtracted from the
right-hand side of the nodal projection, so that :math:`\nabla \cdot U^{n+1} = S^{n+1}`.
