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
  - find real absorption coefficient for CB loadings

![0cb_image_cross-section.png](exports\0cb_image_cross-section.png)
