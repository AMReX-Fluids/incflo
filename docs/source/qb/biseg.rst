.. _Chap:QB:biseg:

Bidisperse Segregation 
====================== 

In this simple test conducted at NETL, a uniform, random mixture 
of nylon 
(:math:`d_p = 3.19` mm, :math:`\rho_p = 1130` kg/m\ :sup:`3` \)
and ceramic 
(:math:`d_p = 4.25` mm, :math:`\rho_p = 2580` kg/m\ :sup:`3` \)
particles. When fluidized, the smaller, lighter nylon particles 
(:math:`U_{mf} \approx 1.1` m/s) segregate out of the mixture to the top, 
with the larger, heavier ceramic particles  (:math:`U_{mf} \approx 1.8` m/s)
remaining at the bottom. The batch segregation test was conducted in 
a small bed with a square cross-section of side length 60.325 mm. 
Similar to previous observations of fluidized segregation of 
bidisperse mixtures [GLMK03]_, the optimal separation was found to 
occur just above the larger of the two  minimum fluidization velocities. 
Below the larger :math:`U_{mf}`, lack of fluidization inhibits particle 
movement and too much above :math:`U_{mf}` vigorous bubbling promotes 
mixture in the bed. 


The optimal batch segregation experiment is simulated with MFiX-Exa 19.08
discretized onto a :math:`8 \times 24 \times 8` mesh. No-slip walls are 
set at vertical domain extents with a mass inflow and pressure outflow at
bottom and top, respectively. The :cpp:`Gidaspow` drag law is applied. 
A defluidization curve of all ceramic particles was traced to find 
:math:`U_{mf} \approx 1.9` m/s, close to the experimental value. The inflow 
velocity is set to `2.0` m/s. The image below shows the simulated segregation
(inset) occurs much more rapidly than observed experimentally. In this case, 
as the nylon particles begin to leave the mixture, the ceramic particles 
defluidize. In the final state, the ceramic particles are essentially static 
with a fluidized layer of nylon particles floating on top, in contrast to the 
experiment, in which the full bed showed signs of fluidization. Investigation 
into the over-segregation in this case is on-going. 
 

.. figure:: figs/netl_biseg_1908_small.png
   :width: 16cm
   :align: center
   :alt: Sim comparison to bidisperse segregation experiment at NETL 

   Comparison of experiment and MFiX-Exa simulaton for rapid segregation
   of a bi-disperse particle mixture.   



