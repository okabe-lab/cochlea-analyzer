function [point2,coeff] = ProjectOntoPCA_Twice(points)
%*****************************************************************************
%This function uses "points", coordinates of points which may represent
%organ of corti. The outputs are "point2", new coordinates of points, and
%"coeff", principal component coefficients.
%*****************************************************************************

%% PCA and noise removal
[coeff,score,~] = pca(points);
[n,~] = size(points);
meanX = mean(points,1);
Xfit = repmat(meanX,n,1) + score(:,1:2)*coeff(:,1:2)';
residuals = points - Xfit;
cut_off = 25;
residual2 = sqrt(residuals(:,1).^2+residuals(:,2).^2+residuals(:,3).^2);
point3 = points(residual2<cut_off,:);

%% PCA and determination of directions of axes
[coeff,~,~] = pca(point3);

meanXYZ = mean(point3);
temp = point3-meanXYZ;
temp_dist = sqrt(temp(:,1).^2+temp(:,2).^2+temp(:,3).^2);
temp_list = temp_dist < median(temp_dist);
meanXYZ2 = mean(point3(temp_list>0,:));
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
point2 = point3*coeff;

f1 = point2(:,2) > median(point2(:,2))+200;
f2 = point2(:,2) < median(point2(:,2))-200;
point2((f1+f2)>0,:)=[];

