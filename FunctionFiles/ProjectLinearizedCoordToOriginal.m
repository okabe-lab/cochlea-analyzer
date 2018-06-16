function xyz2 = ProjectLinearizedCoordToOriginal(xyz1, width, depth, xres, zres ...
    , xres2, zres2, coef, z_mean, rad2, sita2, center, fun_X, coeff)
%*****************************************************************************
%This function uses "xyz1", xyz coordinates in linearized image, "width",
%width of linearized image, "depth", depth of linearized image, "xres" and
%"zres", resolution of original image stack, "xres2" and "zres2" resolution
%of linearized image, "coef" coefficients in equation of fitting curve, 
%"z_mean", mean z coordinates of points, "rad2" and "sita2", angular 
%coordinate and corresponding inclination, "center", pole of polar 
%coordinate system, "fun_X", function handle which transform x coordinate 
%in linearized image into radial coordinate in original image stack. The 
%output is xyz coordinate in original image corresponding to "xyz1".   
%*****************************************************************************

%% Project coordinates in linearlized image into coordinates in original
txyz1 = (xyz1-repmat([0,width+1,depth+1],size(xyz1,1),1))*diag([xres2,xres2,zres2]);
temp_rad = fun_X(txyz1(:,1));

t = temp_rad;
r = [ones(size(t,1),1) t]*coef;
temp_x = r.*cos(t)+center(1);
temp_y = r.*sin(t)+center(2);
temp_z = z_mean;

temp_sita = interp1(rad2,sita2,temp_rad,'linear','extrap');

dxdt = coef(2)*cos(temp_rad)+(coef(2)*temp_rad + coef(1)).*(-sin(temp_rad));
dydt = coef(2)*sin(temp_rad)+(coef(2)*temp_rad + coef(1)).*(cos(temp_rad));
tang_vect = [dxdt dydt]./sqrt(dxdt.^2+ dydt.^2);
rot_vect = [tang_vect zeros(size(xyz1,1),1)];
temp_vect3 = [rot_vect(:,2) -rot_vect(:,1)];
temp_vect = [txyz1(:,2).*temp_vect3 zeros(size(txyz1,1),1)];

txyz2 = zeros(size(xyz1,1),3);
for i = 1:size(xyz1,1)
    temp_vect2 = RodriguesRotation(temp_vect(i,:),rot_vect(i,:),temp_sita(i,:));
    txyz2(i,:) = [temp_x(i,:),temp_y(i,:),temp_z]+temp_vect2;
end

txyz3 = txyz2/coeff;
txyz3(:,3) = txyz3(:,3)+txyz1(:,3);

xyz2 = txyz3*diag([1/xres,1/xres,1/zres]);
