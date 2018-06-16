function xyzOrigin = ConvertLinearToOriginal(xyzLinear,centerLine,inclines,LINEAR_IM_W,LINEAR_IM_D)
%*****************************************************************************
%This function uses "xyzLinear", pixel positions in linearized image,
%"centerLine", point list on center line along hair cells, "inclines",
%vectors indicating proximal directions on approximation planes at each point,
%"LINEAR_IM_W", half size of first dimension in linearized image, and
%"Linear_IM_D", half size of third dimension. The output is "xyzOrigin",
%cartesian coordinates of original image correspoinding to "xyzLinear". 
%*****************************************************************************

centers = interp1(1:size(centerLine,1),centerLine,xyzLinear(:,1),'linear','extrap');
vectors = interp1(1:size(inclines,1),inclines,xyzLinear(:,1),'linear','extrap');
temp = xyzLinear - repmat([0,LINEAR_IM_W+1,LINEAR_IM_D+1],size(xyzLinear,1),1);
xyzOrigin = centers + vectors.*temp(:,2);
xyzOrigin(:,3) = xyzOrigin(:,3) + temp(:,3);