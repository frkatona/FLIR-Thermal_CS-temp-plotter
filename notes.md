for next time:
- read through Holman and note where relevant equations are
 - update powerpoint along the way
- pass over comprehensively for misunderstandings
- append/change any equations I misunderstood previously (e.g., convection that incorporates its previous temperature, power erf function, power with density or no?)
- strip away code until the axes of non-zero graphs make sense again

P slice too high when attenuation is low
Temp higher at bottom even when convection is high

modifiers to control weird distribution (especially conductivity to keep it low)
implementing heat capacity to q addition in T update
flipping the array (or finding where it gets flipped in RK4)

P_slice treated as power vs power density ~line 39 (maybe distribute after P_slice is found)
implement piece-wise nature of thermal conductivity contribution of cb content

changed heat capacity to be in m^3 instead of cm^3

current troubleshooting:
 - convection (top adiabaitic).  Also, top never gets actually hot so its convection isn't accurate
 - conductivity modifier outer = 0.001
- absorption modifier outer = 100
 - Nx = 20
 - flipping q before RK4

additional options: 


q currently in terms of W

current goal:
 - make sure q is just W (Holman p. 145)
 - be sure that the T update uses the volume of space correctly with q


 Book notes to return to:
 [p.31] - radial dimensions for cylinder 
 [p.33] - heat coeff


just trying to get heat to make sense...
current cheat numbers:
(1) conductivity_outer = 1e2
(2) abs_outer = 1e1
(3) P_array *= 1e12

(x) M *= 1e2
(x) RK disabled

fake convection particularly inaccurate with high abs_outer
 - but also, I should be adding q to the convective top I think
 - is there a way to phrase the convective top in terms of the laplacian + q equation(just subtract a term at the end?)

could be helpful to keep track of the total power that has entered and exited the system and include that in the image titles

convection almost fixed, just need to make it perform the heat transfer step too

is my RK4 implementation stepping forward in time or space?

once convection is even remotely passable for transferring heat radially and into the air, start an lmfit attempt where the cheat numbers are parameters
 - can start extracting the surface and depth numbers to fit to...look up how lmfit takes that data (how should I have the simulation output data and how should I extract data from the experimental images in terms of data structures and shapes)
 - can use the first lmfit extraction to get baseline values for things that should be consistent between the sets (I forget which)

power different at different Nx --> must not account for voxel volume appropriately

 QoL:
  - combine preview graphs to 1 window and include in them the total time steps and estimated time