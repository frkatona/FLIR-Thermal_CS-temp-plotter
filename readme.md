# FLIR thermal image-to-temperature-crosssection converter/plotter

- takes folder of FLIR images
- re-maps temperature profiles to a shared scale with custom range and colormap
- plots slice of values across the highest temperature point in each image
- plots log fit of the highest temperature pixel vs time

## Usage

- put all images in a folder
- name images the time in seconds since the start of the experiment (e.g., 0.txt, 10.txt, 60.txt, 300.txt)
- change the file_path variable in the appropriate image_analyzer script and run it

## Data notes

- "low power" is 1.3 W/cm2
- "high power" is 40 W/cm2
- data cut off early when fires started or the fiber optics were becoming damaged

## Relevant values

- thermal conductivity (PDMS) = 0.2 W/mK
- specific gravity (PDMS) = 1.03
- heat capacity, specific (PDMS) = 1 J/gK
- temperature, ambient (air) = 20 C
 
## To-do

- Simulation
- find real absorption coefficient for CB loadings

## Notes

The power \( P(z) \) at a depth \( z \) in a medium due to absorption can be described by the Beer-Lambert law. Given an initial power \( P_0 \) and an absorption coefficient \( \alpha \), the power at a depth \( z \) is given by:

$$
P(z) = P_0 \cdot e^{-\alpha z}
$$

Where:

- \( P(z) \) is the power at depth \( z \).
- \( P_0 \) is the initial power (i.e., power at the surface or \( z = 0 \)).
- \( \alpha \) is the absorption coefficient of the medium.
- \( e \) is the base of the natural logarithm (approximately equal to 2.71828).

![thermal profile - laser - 5 s](exports\upgrade-examples\temp-profile_005s.png)
![thermal profile - laser - 300 s](exports\upgrade-examples\temp-profile-300s.png)
![thermal profile workup - slice](exports\upgrade-examples\slice.png)
![thermal profile workup - log](exports\upgrade-examples\log.png)