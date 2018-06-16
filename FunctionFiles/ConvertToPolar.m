function [radialCoords,angularCoords] = ConvertToPolar(pointList,center)
%*****************************************************************************
%This function uses "pointList", cartesian coordinates of points, "center",
%pole for polar coordinate system. The output is "radialCoords" and
%"angularCoords", polar coordinates correspond to "pointList".
%*****************************************************************************

temp = pointList(:,1:2)-repmat(center,size(pointList,1),1);
radialCoords = sqrt(temp(:,1).^2+temp(:,2).^2);
angularCoords = ComputeAngularCoord(temp);