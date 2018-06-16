function [newpoint,center] = CircleFit_NoiseRemoval2(points,temp_center)
%*****************************************************************************
%This function uses "points", xyz coordinates of points, "temp_center", 
%initial value of center coordinate for circle fitting. The outputs are
%"newpoint", new coordinate of points, "center", center coordinate of
%circle fitted.
%*****************************************************************************

%% Circle fitting
gain = 2:0.2:4;
interval = 0.03:0.01:0.09;
length = zeros(size(gain,2),size(interval,2));
width = zeros(size(gain,2),size(interval,2));

for i = 1:size(gain,2)
    for j = 1:size(interval,2)
        [length(i,j),width(i,j)] = DetectPointsInLine(points,temp_center,gain(i),interval(j));
    end
end

lth = 0.9;
temp = length > lth;
if sum(temp(:)) > 0
    temp2 = (length > lth).*width;
    temp2(temp2==0) = inf;
    [~,idx] = min(temp2(:));
    [I,J] = ind2sub([size(gain,2),size(interval,2)],idx);
else
    [~,idx] = max(length(:));
    [I,J] = ind2sub([size(gain,2),size(interval,2)],idx);
end

[~,~,point2] = DetectPointsInLine(points,temp_center,gain(I),interval(J));
center = ComputeCircleCenter(point2, temp_center, 7);

%% Remove noise
temp = points(:,1:2)-repmat(center,size(points,1),1);
dist = sqrt(temp(:,1).^2+temp(:,2).^2);
mdist = (prctile(dist,20)+prctile(dist,80))/2;
dist2 = dist - ones(size(dist,1),1)*mdist;
f3 = (dist2 < 55).*(dist2 > -55);
newpoint = points(f3>0,:);

