function selectedPoints ...
    = SelectPointsForFit(pointList_PCA,circleCenter,pcCoefficients,imSize,PIXEL_WIDTH,Z_STEP)
%*****************************************************************************
%This function uses "pointList_PCA", points projected in PCA space,
%"circleCenter", pole of polar coordinate system, "pcCoefficients",
%coefficients of principal components, "imSize", size of original image
%stack, "PIXEL_WIDTH", and "Z_STEP". The output is "selectedPoints"
%points selected for further curve fitting.
%*****************************************************************************

[group1_3,center_g1] = CircleFit_NoiseRemoval2(pointList_PCA,circleCenter);
[center_g1, coef_g1] = CircleFitting(group1_3,center_g1);
f_g1 = IsNearBoundaries_Apex(group1_3,center_g1,coef_g1,pcCoefficients,imSize ...
    *diag([PIXEL_WIDTH,PIXEL_WIDTH,Z_STEP]));
group1_4 = group1_3(f_g1==0,:);
selectedPoints = group1_4 / pcCoefficients;