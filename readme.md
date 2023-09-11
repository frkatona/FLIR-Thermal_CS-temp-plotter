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

- find real absorption coefficient for CB loadings
- bugs
  - % transmittance and max np.sum(q) are wild (~ -1e10 and ~ 1e6), probably related to a buggy way of dividing the T and q into the discrete points
    - beam radii below ~0.3 shows 0 q across center as well
  - T distribution never appears isolated to the surface even when q shows a tapering (even with dt = 0.005 and Nx/y/z = 100)