# Microlenses-parameters
Algorithm to analyse parameters of microlenses array based on their hologram images.


## Disclaimer! 
Rest of the project including initial images, files and methods used are not available due to copyright.

## Project goal
The main goal of this project was to develop algorithm to analyse parameters of microlenses based on their hologram images. Maximal automatization of this process was desirable. Project was divided into three parts. At the beginning different numerical propagations algorithms and autofocusing methods was tested to obtain the sharpness possible image of microlenses. The second part was to automatically detect lenses (circles) and get their position coordinates as well as radius of each particular lens. Third and the last part was to connect together two previous tasks, use various phase unwrapping methods and calculate final microlenses parameters.

## General processing path
1. Autofocusing
2. Phase unwrapping
3. Microlenses detection
4. Microlenses parameters calculation

## 1. Autofocusing
Autofocusing For autofocusing angular spectrum (AS) propagation and NORM_L1 method were used. For given propagation range focusing curve was obtained (shown below).

![Autofocusing](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/autofocusing.png)

Focus plane for phase objects like microlenses is obtained in minimum of focus curve. Comparison of initial image and focused one is shown below.

![Hologram plane and image plane](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/hologram%26image_plane.png)

## 2. Phase unwrapping
For phase unwrapping two methods were used: Miguel and TIE. Real calculation time of phase unwrapping was comparable for both methods and due to that MatLab ‘tic toc’ tool was used. Miguel method was more than twice quicker than TIE (9,41875 vs 21,84625 respectively; shown values are average of several runs).

### Before phase unwrapping
![Miguel before phase unwrapping](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/miguel_before.png)
![TIE before phase unwrapping](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/tie_before.png)

### After phase unwrapping
![Miguel after phase unwrapping](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/miguel_after.png)
![TIE](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/tie_after.png)

Both methods provide nearly identical results even on the edges. In this case TIE method can be consider as the better one due to better separation microlenses from the background.

## Microlenses detection (coordinates and radiuses)
In this section only results for TIE method are shown due to very similar results for both methods. After phase unwrapping matix was grayscaled, binarized and incomplete circles from the edges were removed. Results from each step are shown below.

### Grayscaling
![Grayscaling](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/grayscaling.png)

### Binarization
![Binarization](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/binarization.png)

### Image cleaning
![Image cleaning](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/image_cleaning.png)

The next steps were gaussian filtering to avoid sharp edges and edge filtering to obtain only shape of each microlens. Results from each step are shown below.

### Gaussian filtering
![Gaussian filtering](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/gaussian%20filtering.png)

### Edge detection
![Edge detection](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/edge_detection.png)

To find centers coordinates of microlenses for provided range of radiuses values circles were drawn in each bright point from “edge detection” image. Intensity was saved form each iteration in accumulation space. Exemplary images from this process are shown below.

![Radius = 90](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/90.png)
![Radius = 100](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/100.png)
![Radius = 130](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/130.png)

By looking at images above it can be seen that not every circle centers looks the same. This is caused by tilt of initial images plane and it introduces errors.
After that step only areas with the highest intensity in accumulation space were considered as potential circles centers by providing threshold (trial and error method). The next step was voting for centers in which possible values of centers for each circle were averaged. From final data (centers coordinates and radiuses) were plotted (shown below).

### All results
![All results](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/all_results.png)

### Magnification
![Magnification](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/magnification.png)

## 4. Microlenses parameters calculation
Obtained microlenses parameters for both methods and given nominal values are shown in table below.

![TIE results](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/tie.png)

![Miguel results](https://github.com/MateuszPolak/Microlenses-parameters/blob/master/README_imgs/mig.png)

## Conclusions
- Program is NOT fully automatic. It requires form the user to select ROI (region of interest) two times. Also for different input images it will probably require some adjusting.
- Obtained parameters of microlenses are similar to nominal parameters.
- In this case main source of errors was tilt of input image plane. Choosing smaller ROI may help (all calculations was done for maximal number of microlenses on purpose).
- Because these are numerical calculations some errors may occure. Windowing can be reduced by windows smoothening. Also we should avoid unwanted interference due to periodicity od the structure.
- Incomplete Circles and phase artefacts/errors should be removed before further image processing.
