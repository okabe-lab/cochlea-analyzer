function [point2,coeff] = ProjectOntoPCA(points)
%*****************************************************************************
%This function uses "points", coordinates of points which may represent
%organ of corti. The outputs are "point2", new coordinates of points, and
%"coeff", principal component coefficients.
%*****************************************************************************

%% PCA and determination of directions of axes
[coeff,~,~] = pca(points);

meanXYZ = mean(points);
temp = points-meanXYZ;
temp_dist = sqrt(temp(:,1).^2+temp(:,2).^2+temp(:,3).^2);
temp_list = temp_dist < median(temp_dist);
meanXYZ2 = mean(points(temp_list>0,:));
about_center = meanXYZ-meanXYZ2;

temp = cross(coeff(:,1),about_center');
if temp(3,1) < 0
    coeff(:,1) = -coeff(:,1);
end
temp = dot(coeff(:,2),about_center');
if temp < 0
    coeff(:,2) = -coeff(:,2);
end
if coeff(3,3) < 0
    coeff(:,3) = -coeff(:,3);
end
point2 = points*coeff;

f1 = point2(:,2) > median(point2(:,2))+200;
f2 = point2(:,2) < median(point2(:,2))-200;
point2((f1+f2)>0,:)=[];

