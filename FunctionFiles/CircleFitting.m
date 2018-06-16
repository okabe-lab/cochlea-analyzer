function [center,coef,rad,width] = CircleFitting(points,temp_center)
%*****************************************************************************
%This function uses "points", xyz coordinates of points, "temp_center", 
%initial value of center coordinate for circle fitting. The outputs are 
%"center", center of the fitted circle, "coef", radius of fitted circle, 
%"rad", anglular coordinate of points.
%*****************************************************************************

%% Circle fitting
gain = 1.8:0.2:4;
interval = 0.03:0.01:0.10;
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

center = ComputeCircleCenter(point2, temp_center, 9);
width = width(I,J);

%% Calculate radius of fitted circle@
temp = points(:,1:2)-repmat(center,size(points,1),1);
rad = ComputeAngularCoord(temp);
dist = sqrt(temp(:,1).^2+temp(:,2).^2);
number = 7;

int_rad = min(rad):(max(rad)-min(rad))/number:max(rad);
m_dist = zeros(number,1);
m_rad = zeros(number,1);
for i = 1:size(int_rad,2)-1
    f1 = rad>=int_rad(1,i);
    f2 = rad< int_rad(1,i+1);
    f = f1.*f2;
    temp = dist(f>0);
    m_dist(i,1) = (prctile(temp,80) + prctile(temp,10))/2;
    m_rad(i,1) = median(rad(f>0));
end
mdist0 = median(m_dist);
coef = [mdist0 0]';