function [rad,coef_g0,center_g0,points] =  CircleFitting_Alternative(points)
%*****************************************************************************
%This function uses "points", xyz coordinates of points. The outputs are
%"rad", angular coordinate of points, "coef_g0", coefficient of curve
%fitting with angular and radial coorinates of points, "center_g0", center 
%of polar coordinate system, and "points", new coordinates of points.  
%*****************************************************************************

%% Remove noise
zcenter = median(points(:,3));
fz1 = points(:,3) > zcenter+15;
fz2 = points(:,3) < zcenter-15;
fz = fz1 + fz2;
points(fz>0,:) = [];

temp_center = [mean(points(:,1)),mean(points(:,2))+600];
temp = points(:,1:2)-repmat(temp_center,size(points,1),1);
temp_dist = sqrt(temp(:,1).^2+temp(:,2).^2);

dcenter = (prctile(temp_dist,10)+prctile(temp_dist,85))/2;
fd1 = temp_dist > dcenter + 45;
fd2 = temp_dist < dcenter - 45;
fd = fd1 + fd2;
points(fd>0,:) = [];

temp = points(:,1:2)-repmat(temp_center,size(points,1),1);
temp_rad = ComputeAngularCoord(temp);
temp_dist = sqrt(temp(:,1).^2+temp(:,2).^2);

%% Extract representative points
k = convhull(temp_rad,temp_dist);
temp_dist2 = temp_dist(unique(k),:);
temp_th = min(temp_dist2(:,1))+20;
temp_rad2 = temp_rad(unique(k),:);
temp_rad3 = temp_rad2(temp_dist2 < temp_th,:);

inc = 0.075;
num = floor((max(temp_rad3)-min(temp_rad3))/inc);
rp = zeros(num,2);
rep_points = zeros(num,3);
for i = 1:num
    f1 = temp_rad >= (min(temp_rad3) + inc*(i-1));
    if i == num
        f2 = temp_rad <= max(temp_rad3);
    else
        f2 = temp_rad < (min(temp_rad3) + inc*i);
    end
    f = f1.*f2;
    tempr = temp_rad(f>0);
    tempd = temp_dist(f>0);
    
    temp2 = prctile(tempd,5);
    
    num2 = round(prctile(tempd,95)-prctile(tempd,5));
    tempn = zeros(num2,1);
    for j = 1:num2
        tempn(j) = sum(tempd < (temp2 + j));
    end
    
    n1 = round(num2/3)+5;
    n2 = n1*2;
    
    d = diff(tempn(n1:n2,:));
    minf = d == min(d);
    x = sum((n1:n2-1)'.*minf)/sum(minf);
    
    f3 = tempd < (temp2+x);
    
    tempr2 = tempr(f3>0);
    tempd2 = tempd(f3>0);
    
    idx = knnsearch(tempd2,prctile(tempd2,95));
    rp(i,:) = [tempr2(idx) tempd2(idx)];
    
    tp1 = points(f>0,:);
    tp2 = tp1(f3>0,:);
    rep_points(i,:) = tp2(idx,:);
end

%% Curve (Circle) fitting
fun_1 = @(center) CircleFitRSS(rep_points,center);
center_g0 = fminsearch(fun_1,temp_center);

temp = points(:,1:2)-repmat(center_g0,size(points,1),1);
rad = ComputeAngularCoord(temp);

temp2 = rep_points(:,1:2)-repmat(center_g0,size(rep_points,1),1);
temprad = ComputeAngularCoord(temp2);
dist2 = sqrt(temp2(:,1).^2+temp2(:,2).^2);

coef_g0 = [ones(size(temprad,1),1) temprad]\dist2;

