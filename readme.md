# FLIR thermal image-to-temperature-crosssection converter/plotter
 - takes folder of FLIR images
 - re-maps temperature profiles to a shared scale with custom range and colormap
 - plots slice of values across the highest temperature point in each image
 - plots log fit of the highest temperature pixel vs time

## Usage
 - put all images in a folder
 - name images the time in seconds since the start of the experiment (e.g., 0.txt, 10.txt, 60.txt, 300.txt)
 - change the file_path variable in the appropriate image_analyzer script and run it