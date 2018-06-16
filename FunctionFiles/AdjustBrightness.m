function adjustedImage =  AdjustBrightness(image, pointBrightness)
%*****************************************************************************
%This function uses "image", 3D image stack, and "pointBrightness",
%brightness of detected points (estimated hair cells). The output is
%"adjustedImage" 3D image stack with adjustment of brightness.
%*****************************************************************************

backGround = median(image(:));
tempim = (double(image)-backGround)*3000/pointBrightness;
tempim(image < backGround*1.5) = image(image < backGround*1.5);
tempim = tempim + (image>=backGround*1.5)*backGround*1.5;
adjustedImage = uint16(tempim);