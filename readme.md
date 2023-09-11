# FLIR thermal image-to-temperature-crosssection converter/plotter

- takes folder of FLIR images
- re-maps temperature profiles to a shared scale with custom range and colormap
- plots slice of values across the highest temperature point in each image
- plots log fit of the highest temperature pixel vs time

## usage

- put all images in a folder
- name images the time in seconds since the start of the experiment (e.g., 0.txt, 10.txt, 60.txt, 300.txt)
- change the file_path variable in the appropriate image_analyzer script and run it

## notes

[Holman, 2010, p.89]

Heat conduction in finite-difference approximation of a square grid is given by:

$$
T_{m+1, n} + T_{m-1, n} + T_{m, n+1} + T_{m, n-1} - 4T_{m, n} = 0
$$

and heat flux is can be introduced by adding a source term \( q \):

$$
T_{m+1, n} + T_{m-1, n} + T_{m, n+1} + T_{m, n-1} - 4T_{m, n} + \frac{q(\Delta x)^{2}}{k} = 0
$$

[Holman, 2010, p.91]

A solid exposed to a convection boundary condition is computed differently:

$$
-k \Delta y \frac{T_{m, n} - T_{m-1, n}}{\Delta x} - k \frac{\Delta x}{2} \frac{T_{m, n} - T_{m, n+1}}{\Delta y} - k \frac{\Delta x}{2} \frac{T_{m, n} - T_{m, n-1}}{\Delta y} = h \Delta y (T_{m, n} - T_{\infty})
$$

which with x = y reduces the boundary temperature to:

$$
T_{m,n}(\frac{h \Delta x}{k} +2) - \frac{h \Delta x}{k} T_{\infty} - \frac{1}{2}(2T_{m-1, n} + T_{m, n+1} + T_{m, n-1}) = 0
$$

The uniform power distribution across a circular cross section given by a collimated top hat distribution yields:

$$
P_A = \frac{P_z}{\pi r^2}
$$

where \( P \) is the power at a given volume element, \( P_{in} \) is the total power input, and \( r \) is the radius of the beam.

The power \( P(z) \) at a depth \( z \) described by the Beer-Lambert law. Given an initial power \( P_0 \) and an absorption coefficient \( \alpha \), the power at a depth \( z \) is given by:

$$
P(z) = P_0 \cdot e^{-\alpha z}
$$



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

- bugs
  - beam mask not being appropriately applied to the T updates
  - possibly bugged how I'm dividing the T and q into the discrete points (weird when x/y are even vs odd at certain beam radii for finding centerline, even floored)
    - % transmittance and max np.sum(q) are bugged (~ -1e10 and ~ 1e6)
  - T distribution never appears isolated to the surface even when q shows a tapering (even with dt = 0.005 and Nx/y/z = 100)
- find real absorption coefficient for CB loadings