function [width,temp_center] = SimpleCircleFitEvaluation(point4,r)
%*****************************************************************************
%This function uses "point4", xy coordinates of points (after PCA), "r", 
%radius of circle. The outputs are "width", expanse of points in radial
%direction, and "temp_center", coordinates of center of circle.
%*****************************************************************************

f1 = point4(:,1) > prctile(point4(:,1),33);
f2 = point4(:,1) < prctile(point4(:,1),67);
point4_2 = point4((f1.*f2)>0,:);
temp_center = [mean(point4_2(:,1)) mean(point4_2(:,2))+r];

temp = point4(:,1:2)-repmat(temp_center,size(point4,1),1);
dist = sqrt(temp(:,1).^2+temp(:,2).^2);

f3 = f1==0;
f4 = f2==0;

width = std([prctile(dist(f3>0),20) prctile(dist(f4>0),20) prctile(dist(f1.*f2>0),20)]);
