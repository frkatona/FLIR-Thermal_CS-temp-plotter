# FLIR thermal image-to-temperature-crosssection converter/plotter

- takes folder of FLIR images
- re-maps temperature profiles to a shared scale with custom range and colormap
- plots slice of values across the highest temperature point in each image
- plots log fit of the highest temperature pixel vs time

## using the FLIR image analyzer

- put all images in a folder
- name images the time in seconds since the start of the experiment (e.g., 0.txt, 10.txt, 60.txt, 300.txt)
- change the file_path variable in the appropriate image_analyzer script and run it

## unsteady-state numerical considerations (*J.P. Holman*, 2010)

### Conduction Through Bulk

In the finite differences method for **unsteady-state** heat flow systems, temperature ($T$) can be found at a given node $(m, n)$ and step ($p$) after a time increment $\Delta \tau$ for spatial increments $\Delta x$ and $\Delta y$ in a material with thermal diffusivity $\alpha$ through:

$$
\frac{T_{m+1, n}^p + T_{m-1, n}^p - 2T_{m, n}^p}{\Delta x^2} + \frac{T_{m, n+1}^p + T_{m, n-1}^p - 2T_{m, n}^p}{\Delta y^2} = \frac{1}{\alpha}\frac{T_{m, n}^{p+1} - T_{m, n}^p}{\Delta \tau} 
$$

And where $\Delta x = \Delta y$, the equation for $T_{m,n}^{p+1}$ becomes:

$$
\mathbf{(1)}\space{}T_{m,n}^{p+1} = \frac{\alpha \Delta \tau}{(\Delta x)^2} (T_{m+1,n}^p + T_{m-1,n}^p + T_{m,n+1}^p + T_{m,n-1}^p) + \left[1-\frac{4\alpha \Delta \tau}{(\Delta x)^2}\right] T_{m,n}^p 
$$

### Heat Source
In the case of constant heat flux on a semi-finite solid, boundary conditions can be described as:

$$
T(x, 0) = T_i
$$

$$
\frac{q_0}{A} = -k \frac{\partial T}{\partial x} \bigg|_{x=0}
$$

and the solution in this case is given by [Holman p. 145]:

$$
T - T_i = \frac{2q_0\sqrt{\alpha \tau / \pi}}{kA} exp(\frac{-x^2}{4\alpha \tau})- \frac{q_0x}{kA} (1-erf\left[\frac{x}{2\sqrt{\alpha \tau}}\right])
$$

[CONSIDER if the erf term applies to the numerical solution]

[CONSIDER how the volumetric nature of the space steps must be taken into account for possible theta values (should dx be chosen based on $\frac{2\pi}{n}$?) and finding the number of nodes that Q will be distributed across]

[CONSIDER if boundary conditions should be additive for the convective and heat source overlap region]

[CONSIDER if the volumetric heat application in eq2 is correct]

[CONSIDER if y integrations should be from y->(y+dy) or (y-1/2 dy)->(y+1/2 dy)]

For a constant power source, a power term can be added to the node $(m, n)$ in $\mathbf{eq1}$.  For the volume of the material which overlaps with the power source, $T_{m,n}^{p+1}$ can be described by:

$$
\mathbf{(2)}\space{}T_{m,n}^{p+1} = \frac{\alpha \Delta \tau}{(\Delta x)^2} (T_{m+1,n}^p + T_{m-1,n}^p + T_{m,n+1}^p + T_{m,n-1}^p) + \left[1-\frac{4\alpha \Delta \tau}{(\Delta x)^2}\right] T_{m,n}^p + \frac{q}{c \rho} \Delta \tau
$$

where $c$ is the specific heat capacity and $\rho$ is the density of the material (combining to describe the volumetric heat capacity) and $q$ is the volumetric power which, in the simplest case, is given by:

$$
q = \frac{Q}{\Delta x \cdot \Delta y \cdot \Delta z}
$$

where, for our purposes, $\Delta x$ is the radial coordinate, $\Delta y$ is axial (i.e., height), and $\Delta z$ is tangent to the curve (note: it's possible an azimuthal angle $\theta$ would be more appropriate--which would yeild cylindrical-shell volume elements for the purposes of distributing heat in *Consideration 2* below--but for now, the volume elements are cubic). However, the power source is not perfectly uniformly distributed across the given volume of the material.  The following three Considerations must be made:

1) depth distribution 
2) radial distribution
3) transmittance

(1) Because the heat source is modeling the absorption of light, the depth distribution of the laser power (i.e., vertical cross-sections) can be described by the Beer-Lambert law.  This exponential decay of power across the depth of the material, $Q(z)$ is given by:

$$
Q_y = Q_0 \cdot e^{-\alpha y}
$$

where $Q_0$ is the power at the surface and $\alpha_{abs}$ is the absorption coefficient (not to be confused with the thermal diffusivity, $\alpha$).

Finding the $q$ for a given voxel requires choosing a $Q_y$ value to represent the voxel which will become less accurate with increasing voxel dimensions.  To mitigate this, the $Q_y$ can be found by integrating the power distribution across the voxel height, shown here:

$$
Q_y = \int_{y}^{y + \Delta y} Q_0 \cdot e^{-\alpha y} dy
$$

(2) Due to the top-hat distribution of the beam in this particular case, the radial/azimuthal distribution (i.e., horizontal cross-sections) is indeed roughly homogeneous/isotropic and therefore independent of x and z (within the bounds of the beam). So, at a given depth, $y$, the area power $Q_{\Delta x \Delta z}$ can be described by its fraction of the circular cross section across the cylinder:

$$
Q_{\Delta x\Delta z} = Q_y\frac{\Delta x \cdot \Delta z}{\pi r^2} 
$$

which makes the power density, $q$, for a given voxel:

$$
q_{xyz} = \left[\int_{y}^{y + \Delta y} Q_0 \cdot e^{-\alpha y} dy \right] \cdot \left[\frac{\Delta x \cdot \Delta z}{\pi r^2}\right]
$$

(3) It is also not necessarily the case that all of the power is converted to heat within the material.  For sufficiently low $\alpha _{abs}$, non-negligible power will be transmitted through the length of the material.  Given the above formulations for q, this should not impact the simulation accuracy, but for any other purposes, the transmittance through the entire height, $h$, of the material can be found by:

$$
\%T = \frac{\int_{h}^{\infty} Q \cdot e^{-\alpha y} dy}{Q} \cdot 100\%
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

which neglects the heat capacity of the element of the wall at the boundary, but that should be small when a large number of $\Delta x$ nodes are used such that the surface area is small relative to the bulk.  Heat capacity can be accounted for, however (when $\Delta x = \Delta y$) with:

$$
\mathbf{(3)}\space{}T_{m,n}^{p+1} = \frac{\alpha \Delta \tau}{(\Delta x)^2} \left[ 2\frac{h\Delta x}{k} T_{\infty} + 2T_{m-1,n}^p + T_{m,n+1}^p + T_{m,n-1}^p + \left[\frac{(\Delta x)^2}{\alpha \Delta \tau} - 2\frac{h \Delta x}{k} - 4\right]T_{m,n}^p \right]
$$

### Satisfying the CFL Condition

For all such systems, convergence relies on a specific relationship between time and space increments and the velocity magnitude ($\alpha$ in this case), known as the Courant-Friedrichs-Lewy (CFL) condition.  For 2D cases such as those presented in equations 1 and 2, the following relationship assures compliance with the CFL condition:

$$
M = \frac{(\Delta x)^2}{\alpha \Delta \tau} \geq 4
$$

Under this condition, $T$ at node $(m, n)$ after $\Delta \tau$ is simply the average of the four surrounding nodes.  At the convective boundary in equation 3, the CFL condition is given by:

$$
\frac{(\Delta x)^2}{\alpha \Delta \tau} >= 2 (\frac{h \Delta x}{k} + 1)
$$

### Combining all Three Cases

And so with each 

## implementing DSC data for cure profiles

asdf

## fitting penetration and T profiles to model

asdf

## physical parameter notes

- "low power" is 1.3 W/cm2
- "high power" is 40 W/cm2
- thermal conductivity (PDMS) = 0.2 W/mK
- specific gravity (PDMS) = 1.03
- heat capacity, specific (PDMS) = 1 J/gK
- temperature, ambient (air) = 20 C
- data cut off early when fires started or the fiber optics were becoming damaged

## example images

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
  - T distribution never appears to decay into the depth though q shows a tapering (even with $\Delta \tau$ = 0.005 and Nx/y/z = 100 and at higher t values)
- find real absorption coefficient for CB loadings

## misc
### Extending to 3D with Cylindrical Symmetry (ignore for now)
The cylindrical power source of the laser beam will require some 3D considerations at least insofar as the proper distribution of the laser power as a heat source is concerned, as well as any 3D visualization of the temperature distribution.  The heat equation in cylindrical coordinates is given by:

$$
\frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial T}{\partial r} \right) + \frac{1}{r^2} \frac{\partial^2 T}{\partial \theta^2} + \frac{\partial^2 T}{\partial z^2} = \frac{1}{\alpha} \frac{\partial T}{\partial t}
$$

Due to the symmetry of the system, the temperature is independent of the azimuthal angle, $\theta$, and so the heat equation reduces to:

$$
\frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial T}{\partial r} \right) + \frac{\partial^2 T}{\partial z^2} = \frac{1}{\alpha} \frac{\partial T}{\partial t}
$$

Which is described by the 2D numerical solutions above. 


- [p. 95] Large number of nodes unnecessary due to inherently large uncertanties in h
- [p. 99] Guass-Seidel iterative method convergence
- [p. 102] Error of the finite-difference approximation to ∂T/∂x is of the orderof (\Delta x/L)^2 where L is some characteristic body dimension (but also floating point errors increase with the number of iterations)
