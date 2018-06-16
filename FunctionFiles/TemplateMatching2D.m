function correlationMat = TemplateMatching2D(targetIm, templateIm)
%*****************************************************************************
%This function uses "targetIm", 2D image, and "templateIm", 2D template
%image. The output is matrix of correlation coefficients.
%*****************************************************************************

pad1 = floor(size(templateIm,1)/2);
pad2 = floor(size(templateIm,2)/2);
imSize = size(targetIm);
correlationMat = normxcorr2(templateIm,targetIm);
correlationMat = correlationMat(pad1+1:pad1+imSize(1),pad2+1:pad2+imSize(2));
