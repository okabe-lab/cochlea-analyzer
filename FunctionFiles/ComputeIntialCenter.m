function initialCenter = ComputeIntialCenter(pointList)
%*****************************************************************************
%This function uses "pointList", cartesian coordinates list. The outputs is 
%"initialCenter", initial value of center coordinate for circle fitting.
%*****************************************************************************

% Estimate center coordinate of circle with limited range of radius
radius = (360:20:800)';
radiusNo = size(radius,1);
radialCoordRange = zeros(radiusNo,1);
tempCenter = zeros(radiusNo,2);
for i = 1:radiusNo
    [radialCoordRange(i,1),tempCenter(i,:)] ...
        = SimpleCircleFitEvaluation(pointList,radius(i,1));
end
[~,idx] = min(radialCoordRange);
initialCenter = tempCenter(idx,:);