# Micronuclei lamina analysis (Mammel et al. 2021)

## MANUSCRIPT

Chromosome length and gene density contribute to micronuclear membrane stability
Anna E Mammel, Heather Z Huang, Amanda L Gunn, Emma Choo, Emily M Hatch
BioRxiv Preprint: https://www.biorxiv.org/content/10.1101/2021.05.12.443914v1

## SUMMARY

This repository contains MATLAB functions designed to segment the lamin network of micronuclei and to identify gaps in the lamin mesh. It has been optimized for image datasets acquired on a Leica SP8 instrument equipped with a STED module. The image stack has 5 different channels. The main function `LaminDetector` segments the lamin network of micronuclei using the AROS algorithm developped by Mark Kittisopikul and colleagues to identify gaps in the mesh from 3D image datasets acquired with super res STED technique.

1. After image import, a montage of all markers is created to help identify the micronucleus of interest.
2. A new interactive figure pops up for cropping the image stack down to the region of interest.
3. In case the micronucleus of interest cannot be properly segmented (e.g. it touches the main nucleus), a new interactive figure opens up to allow the user to circle the object of interest.
4. The lamin network is then segmented in 3D using the AROS algorithm, binarized, split into a top and a bottom hemispheres and holes properties within those hemispheres are extracted.
5. A final figure showing the relevant markers, the intensity-weighted segmented lamin mesh, and the identified holes in the top and bottom hemispheres is generated for review.

## REQUIREMENTS

MATLAB R2019a or later

    Image Processing Toolbox
    Statistics and Machine Learning Toolbox
    Parallel Computing Toolbox

`bfmatlab` package, the MATLAB bioformat image converter that can be downloaded at:
https://docs.openmicroscopy.org/bio-formats/6.3.1/users/matlab/

`AdaptiveResolutionOrientationSpace` that can be found at:
https://github.com/mkitti/AdaptiveResolutionOrientationSpace

`fmeasure` , a MATLAB function to measure image focus:
https://www.mathworks.com/matlabcentral/fileexchange/27314-focus-measure

`ellipsoid_fit` , a MATLAB function to fit an ellipsoid to a point cloud:
https://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit


