Creating the MAC velocities
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To create the normal velocities on faces, we first extrapolate from the cell centers on each side using the
slopes as computed earlier, and upwind the face value to define  :math:`U^{pred}`

.  To compute the x-velocity on the x-faces of regular (ie not cut) cells, we call

   .. code:: shell

            AMREX_CUDA_HOST_DEVICE_FOR_3D(ubx, i, j, k,
             {
                 // X-faces
                 Real upls     = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                 Real umns     = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                 if ( umns < 0.0 && upls > 0.0 ) {
                    umac_fab(i,j,k) = 0.0;
                 } else {
                    Real avg = 0.5 * ( upls + umns );
                    if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                    } else if (avg >= 0)             { umac_fab(i,j,k) = umns;
                    } else                           { umac_fab(i,j,k) = upls;
                    }
                 }
             });

For cut cells we test on whether the area fraction is non-zero: 
       
   .. code:: shell

             AMREX_CUDA_HOST_DEVICE_FOR_3D(ubx, i, j, k,
             {
                 // X-faces
                 if (ax_fab(i,j,k) > 0.0)
                 {
                    Real upls     = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                    Real umns     = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                    if ( umns < 0.0 && upls > 0.0 ) {
                       umac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( upls + umns );
                       if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { umac_fab(i,j,k) = umns;
                       } else                           { umac_fab(i,j,k) = upls;
                       }
                    }
                 } else {
                       umac_fab(i,j,k) = huge_vel;
                 }
             });

We then perform a MAC projection on the face-centered velocities to enforce that they satisfy 

.. math:: \nabla \cdot (U^{MAC})  = 0

We do this by solving 

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi^{MAC} = \nabla \cdot \left(U^{pred} \right)

then defining

.. math:: U^{MAC} = U^{pred} - \frac{1}{\rho} \nabla \phi^{MAC}
