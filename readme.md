# FLIR thermal image-to-temperature-crosssection converter/plotter

- takes folder of FLIR images
- re-maps temperature profiles to a shared scale with custom range and colormap
- plots slice of values across the highest temperature point in each image
- plots log fit of the highest temperature pixel vs time

## usage

- put all images in a folder
- name images the time in seconds since the start of the experiment (e.g., 0.txt, 10.txt, 60.txt, 300.txt)
- change the file_path variable in the appropriate image_analyzer script and run it

## numerical solution considerations (*J.P. Holman*, 2010)

### Steady State

In steady state systems, **Thermal conduction** in finite-difference approximation of a square grid is given by [Holman p.89]:

$$
T_{m+1, n} + T_{m-1, n} + T_{m, n+1} + T_{m, n-1} - 4T_{m, n} = 0
$$

and a **heat source** is can be introduced by adding a source term \( q \):

$$
T_{m+1, n} + T_{m-1, n} + T_{m, n+1} + T_{m, n-1} - 4T_{m, n} + \frac{q(\Delta x)^{2}}{k} = 0
$$

A solid exposed to a **convection boundary** condition is computed differently [Holman, p.91]:

$$
-k \Delta y \frac{T_{m, n} - T_{m-1, n}}{\Delta x} - k \frac{\Delta x}{2} \frac{T_{m, n} - T_{m, n+1}}{\Delta y} - k \frac{\Delta x}{2} \frac{T_{m, n} - T_{m, n-1}}{\Delta y} = h \Delta y (T_{m, n} - T_{\infty})
$$

which with **x = y reduces** the boundary temperature to:

$$
T_{m,n}(\frac{h \Delta x}{k} +2) - \frac{h \Delta x}{k} T_{\infty} - \frac{1}{2}(2T_{m-1, n} + T_{m, n+1} + T_{m, n-1}) = 0
$$

### Unsteady State

In an **unsteady-state** system:

$$
\frac{T_{m+1, n}^p + T_{m-1, n}^p - 2T_{m, n}^p}{\Delta x^2} + \frac{T_{m, n+1}^p + T_{m, n-1}^p - 2T_{m, n}^p}{\Delta y^2} = \frac{1}{\alpha}\frac{T_{m, n}^{p+1} - T_{m, n}^p}{\Delta t}
$$

and for such 2-D systems, in situations where time and distance increments are chosen as:

$$
M = \frac{(\Delta x)^2}{\alpha \Delta \tau} = 4
$$

then T at node (m, n) after dt is simply the average of the four surrounding nodes.  

Known as the CFL condition, this also reveals the values of dt that can feasibly be chosen for a given dx.  When M is low, the equation allows heat to flow from a cold node to a hot node, and so accuracy relies on the condition:

$$
M \geq 4
$$

### Convection Boundary

Above relations do not apply at convection boundaries and must be handled separately.  In the case of a flat wall in 2-D, the finite-difference approximation is given by [Holman p. 170]:

$$
-k\frac{\Delta y}{\Delta x} (T_{m+1} - T_{m}) = h \Delta y (T_{m+1} - T_{\infty})
$$

or rearranged to

$$
T_{m+1} = \frac{T_m + (h \Delta x / k) T_{\infty}}{1 + h \Delta x / k}
$$

which neglects the heat capacity of the element of the wall at the boundary, but that should be small when a large number of dx nodes are used such that the surface area is small relative to the bulk.  Heat capacity can be accounted for, however (when dx = dy) with:

$$
T_{m,n}^{p+1} = \frac{\alpha \Delta \tau}{(\Delta x)^2} \left[ 2\frac{h\Delta x}{k} T_{\infty} + 2T_{m-1,n}^p + T_{m,n+1}^p + T_{m,n-1}^p + \left[\frac{(\Delta x)^2}{\alpha \Delta \tau} - 2\frac{h \Delta x}{k} - 4\right]T_{m,n}^p \right]
$$

and to ensure convergence of the numerical solution, selections of dx and dt must satisfy the convective 2D CFL condition:

$$
\frac{(\Delta x)^2}{\alpha \Delta \tau} >= 2 (\frac{h \Delta x}{k} + 1)
$$
### Power terms
### Cylindrical Symmetry System
The uniform power distribution across a circular cross section given by a collimated top hat distribution yields:

$$
P_A = \frac{P_z}{\pi r^2}
$$

where \( P \) is the power at a given volume element, \( P_{in} \) is the total power input, and \( r \) is the radius of the beam.

The power \( P(z) \) at a depth \( z \) described by the Beer-Lambert law. Given an initial power \( P_0 \) and an absorption coefficient \( \alpha \), the power at a depth \( z \) is given by:

$$
P(z) = P_0 \cdot e^{-\alpha z}
$$


- [p. 95] Large number of nodes unnecessary due to inherently large uncertanties in h
- [p. 99] Guass-Seidel iterative method convergence
- [p. 102] Error of the finite-difference approximation to ∂T/∂x is of the orderof (dx/L)^2 where L is some characteristic body dimension (but also floating point errors increase with the number of iterations)

### physical parameter notes

- "low power" is 1.3 W/cm2
- "high power" is 40 W/cm2
- thermal conductivity (PDMS) = 0.2 W/mK
- specific gravity (PDMS) = 1.03
- heat capacity, specific (PDMS) = 1 J/gK
- temperature, ambient (air) = 20 C
- data cut off early when fires started or the fiber optics were becoming damaged

### thermal profile - laser - 5 s

![thermal profile - laser - 5 s](exports\upgrade-examples\temp-profile_005s.png)

### thermal profile - laser - 300 s

![thermal profile - laser - 300 s](exports\upgrade-examples\temp-profile-300s.png)

### thermal profile workup - slice

![thermal profile workup - slice](exports\upgrade-examples\slice.png)

### thermal profile workup - log

![thermal profile workup - log](exports\upgrade-examples\log.png)

## to-do

- refactor simulation for cylindrical symmetry
- investigate bugs
  - beam mask not being appropriately applied to the T updates
  - npz and print temperatures don't reflect the temperatures  in the graph
  - transmittance shows % <<0 and >>100 (problem with how I'm discretizing power?)
  - T distribution never appears to decay into the depth though q shows a tapering (even with dt = 0.005 and Nx/y/z = 100 and at higher t values)
- find real absorption coefficient for CB loadings
