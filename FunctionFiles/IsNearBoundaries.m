function [tf,fit_curve_g3] = IsNearBoundaries(points, center, coef1, coef2, Isize)
%*****************************************************************************
%This function uses "points", xyz coordinates of points, "center", pole of 
%polar coordinate, "coef1", coefficient of fitted curve with radial and 
%polar coodinates, "coef2", principal component coefficients, and "Isize", 
%size of image stack. The outputs are "tf", result of selection of points,
%"fit_curve_g3", xyz coordinates of fitting curve.
%*****************************************************************************

%% Remove points near boundaries of image stack
temp = points(:,1:2)-repmat(center,size(points,1),1);
rad = ComputeAngularCoord(temp);

coef0 = coef1;
coef3 = coef1;

coef1(1) = coef1(1)+25;
r = [ones(size(rad,1),1) rad]*coef1;
xx1 = r.*cos(rad)+center(1);
yy1 = r.*sin(rad)+center(2);
fit_curve_g1 = [xx1 yy1 mean(points(:,3))*ones(size(xx1,1),1)];

fit_curve_g1 = fit_curve_g1 / coef2;
tf = (fit_curve_g1(:,1) < 0)+(fit_curve_g1(:,1) > Isize(1))+(fit_curve_g1(:,2) < 0) ...
    +(fit_curve_g1(:,2) > Isize(2));

coef0(1) = coef0(1)-40;
r = [ones(size(rad,1),1) rad]*coef0;
xx1 = r.*cos(rad)+center(1);
yy1 = r.*sin(rad)+center(2);
fit_curve_g2 = [xx1 yy1 mean(points(:,3))*ones(size(xx1,1),1)];
fit_curve_g2 = fit_curve_g2 / coef2;

tf = tf + (fit_curve_g2(:,1) < 0)+(fit_curve_g2(:,1) > Isize(1))+(fit_curve_g2(:,2) < 0) ...
    +(fit_curve_g2(:,2) > Isize(2));
tf = tf + (fit_curve_g2(:,3) < 0)+(fit_curve_g2(:,3) > Isize(3));

tf = tf + (abs(points(:,3)-mean(points(:,3)))>20);
tf = tf>0;

r = [ones(size(rad,1),1) rad]*coef3;
xx1 = r.*cos(rad)+center(1);
yy1 = r.*sin(rad)+center(2);
fit_curve_g3 = [xx1 yy1 mean(points(:,3))*ones(size(xx1,1),1)];
fit_curve_g3 = fit_curve_g3 / coef2;

if tf(round(numel(tf)/2)) == 1
    tf2 = bwlabel(tf);
    tf3 = tf2==tf2(round(numel(tf)/2));
    tf(tf3>0) = 0;
end

