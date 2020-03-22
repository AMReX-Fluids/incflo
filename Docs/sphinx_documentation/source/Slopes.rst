Computing slopes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Slopes are computed one direction at a time for each scalar and each
velocity component.

We use the second order Monotonized Central (MC)
limiter (van Leer, 1977). The scheme is described below for the u-velocity.

The limiter computes the slope at cell "i" by combining the left, central
and right u-variation "du":

.. code:: shell

     du_l = u(i) - u(i-1)               = left variation
     du_c = 0.5 * ( u(i+1) - u(i-1) )   = central (umlimited) variation
     du_r = u(i+1) - u(i)               = right variation

Finally, the u-variation at cell "i" is given by :

.. code:: shell

     du(i) = sign(du_c) * min(2|du_l|, |du_c|, 2|du_r|)) if du_l*du_r > 0
     du(i) = 0                                           otherwise

The above procedure is applied direction by direction.

BOUNDARY CONDITIONS
When periodic or Neumann's BCs are imposed, the scheme can be applied
without any change since the ghost cells at the boundary are filled
by either periodicity or by extrapolation.
For Dirichlet's BCs in the transversal direction, the scheme can again
be applied as is since the velocity is known at the first ghost cell
out of the domain.
However, for Dirichlet's BCs in the longitudinal direction, the velocity
is not known outside the domain since the BC is applied directly at the first
valid node which lies on the boundary itself. Therefore, the scheme must be
arranged as follows to use ONLY values from inside the domain.
For a left boundary (i=0), the u-variations are:

.. code:: shell

     du_l = 0                             Dont use values on the left
     du_c = -1.5*u(0) + 2*u(1) -0.5*u(2)  2nd order right-biased
     du_r = u(1) - u(0)                   Right variation
