function isInIm = IsWithinIm(pixelList,imSize)
%*****************************************************************************
%This functin uses "pixelList", 3D pixel positions, "imSize", image size.
%The output is "isInIm", logical matrix indicating within image 
%boundaries (1) or not (0).
%*****************************************************************************

inIm1 = (pixelList(:,1) >= 1) .* (pixelList(:,1) <= imSize(1));
inIm2 = (pixelList(:,2) >= 1) .* (pixelList(:,2) <= imSize(2));
inIm3 = (pixelList(:,3) >= 1) .* (pixelList(:,3) <= imSize(3));
isInIm = inIm1.*inIm2.*inIm3;