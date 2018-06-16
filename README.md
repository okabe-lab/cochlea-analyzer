
Matlab scripts for analyzing image stacks of cleared cochlea
_____________________________________________________________________________________________

This program needs:
  - Matlab (R2017b or newer version)
  - Image Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Neural Network Toolbox

Data size:
  Total         514 MB
    - Program    77 MB
    - Test data 437 MB
_____________________________________________________________________________________________

For testing the program, please run the following scripts from top to bottom:
  - main0.m (unzip test data, 1 min.)
  - main1.m (make linearized image, 15 min.)
  - main2.m (detect inner hair cells, 5 min.)
  - main3.m (detect outer hair cells, 3 min.).
  - main4.m (analyze spatial distribution of outer hair cells, <1 min.) 
  - main5.m (make cellular cartography, <1 min.)
  - main6.m (simulation of clustered cell loss, 10 min.)

It takes ~30 minutes for running all the scripts on PC with:
   Windows 10 Home,
   Inter Core i7-6700 CPU @ 3.40GHz,
   16.0GB RAM.

Following result files will be saved in "..\TestData\Result":
  - "linearizedIm1.tif", "linearizedIm2.tif"     (Linearized images of organ of Corti)
  - "innerHairCells.tif", "outerHairCells.tif"   (Images indicating detected cells)
  - "innerHairCells.xlsx", "outerHairCells.xlsx" (Coordinates of detected cells)
  - "standardized_OuterHairCells.tif"            (Standardized image of cellular distribution) 
_____________________________________________________________________________________________

June 16, 2018

Tadatsune Iida
Department of Cellular Neurobiology
Graduate School of Medicine, The University of Tokyo
Email: tadiida@m.u-tokyo.ac.jp

