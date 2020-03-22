Fluid Variables
===============

   +-----------------------+--------------------------------------------------+
   | Variable              | Definition                                       |
   +=======================+==================================================+
   | :math:`\rho_g`        | Fluid density                                    |
   +-----------------------+--------------------------------------------------+
   | :math:`\varepsilon_g` | Volume fraction of fluid (= 1 if no particles)   |
   +-----------------------+--------------------------------------------------+
   | :math:`U_g`           | Fluid velocity                                   |
   +-----------------------+--------------------------------------------------+
   | :math:`\tau`          | Viscous stress tensor                            |
   +-----------------------+--------------------------------------------------+
   | :math:`g`             | Gravitational acceleration                       |
   +-----------------------+--------------------------------------------------+

Fluid Equations
===============

Conservation of fluid mass:

.. math:: \frac{\partial (\varepsilon_g \rho_g)}{\partial t} + \nabla \cdot (\varepsilon_g \rho_g U_g)  = 0

Unlike two-fluid modeling, :math:`\varepsilon_g = 1 - \varepsilon_s`, is not a solution 
variable. Rather, :math:`\varepsilon_s` is evaluated from the particle field through 
volume filtering, 

.. math:: 
   (1 - \varepsilon_g) A (\boldsymbol{x},t) \approx 
   \sum_{i=1}^{N_p} A_i(\boldsymbol{X}_i,t) {\mathcal G} 
   (\left| \boldsymbol{x} - \boldsymbol{X}_i \right|) {\mathcal V}_i 

where :math:`A_i` is a general particle property, :math:`{\mathcal V}_i` is the particle 
volume and :math:`\mathcal G` is a transfer kernel with compact support--here linear hat. 
Setting :math:`A_i = 1` gives the local void fraction.  


Assuming the fluid phase is incompressible, :math:`\frac{D \rho_g}{Dt} = 0`, the conservation of fluid mass is equivalent to the conservation of fluid volume:

.. math:: \frac{\partial \varepsilon_g}{\partial t} + \nabla \cdot (\varepsilon_g  U_g)  = 0

The conservation of fluid momentum is:

.. math:: \frac{ \partial (\varepsilon_g \rho_g U_g)}{\partial t} 
   + \nabla \cdot (\varepsilon_g \rho_g U_g U_g) + \varepsilon_g \nabla p_g
   = \nabla \cdot \tau + M_{sg} + \varepsilon_g \rho_g g

where :math:`M_{sg} = - M_{gs}` is the generalized interfacial momentum transfer from 
the solid particles to the fluid-phase. Like :math:`\varepsilon_s`, 
:math:`M_{gs}` is determined from the L-E transfer kernel by setting :math:`A_i = F_{gi}`, 
where :math:`F_{gi}` is the force due to the fluid-phase on the ith particle. Following 
MFiX classic (and many other CFD-DEM codes designed for high density ratio gas-solids 
flows), only buoyancy (pressure gradient) and steady drag are considered: 

.. math::
   F_{gi} = - \mathcal{V}_i \nabla p_g 
   - \frac{1}{2} C_D \rho_g \boldsymbol{V}_{ig} \left|\boldsymbol{V}_{ig}\right| A_i^{(proj)}

where :math:`\boldsymbol{V}_{ig} = \boldsymbol{V}_i - \boldsymbol{U}_g ( \boldsymbol{X}_i )` 
is the velocity of ith particle relative to the fluid-phase (at the particle position 
:math:`\boldsymbol{X}_i`). :math:`F_{gi}` is closed by the specification of a drag 
coefficient, :math:`C_D`. Currently, MFiX-Exa includes Wen-Yu, Gidaspow and BVK2 drag laws.

.. wdf todo 
   gidaspow form: C_D^{Gidaspow} = \chi C_D^{(Wen-Yu)} + \left(1 - \chi \right)  C_D^{(Ergun)}
   \chi = \frac{\arctan 150 \left( \varepsilon_g - 0.8 \right)}{\pi} + \frac{1}{2}
   wen-yu:        C_D^{(Wen-Yu)} = \max \left[\frac{24}{Re_i}\left(1 + 0.15 Re_i^{0.687}\right),\ 0.44 \right] \left( 1 - \varepsilon_g \right)^{-1.65}
   ergun: C_D^{(Ergun)} = \frac{200 \left(1 - \varepsilon_g\right)}{Re_i} + \frac{7}{3} 
   BVK2: C_D^{BVK2} = ...
   particle reynolds number: Re_i = \frac{\rho_g \left( 1 - \varepsilon_g \right) d_i \left| \boldsymbol{V}_{ig} \right| }{\mu_g}


In chemical engineering literature, it is common to lump all drag-related terms of 
:math:`M_{gs}` into :math:`\beta`. With this simplification and some re-arrangement, 
the fluid momentum takes the more convenient form: 

.. math:: \frac{ \partial (\varepsilon_g \rho_g U_g)}{\partial t} 
   + \nabla \cdot (\varepsilon_g \rho_g U_g U_g) + \varepsilon_g \nabla p_g 
   = \nabla \cdot \tau + \sum_p \beta_p (V_p - U_g) + \varepsilon_g \rho_g g



