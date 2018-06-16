function [centerline1, inclines] = DrawLinesAlongCorti(o_point, spoint, PIXEL_WIDTH, Z_STEP)
%*****************************************************************************
%This function uses "o_point", xyz coordinates of all points, "spoint",
%points selected for curve fitting, and "PIXEL_WIDTH" and "Z_STEP", image
%resolution. The outputs are "centerline1", fitting curve, and
%"centerline2", fitting curve slightly shifted to proximal.
%*****************************************************************************

%% PCA and determination of directions of axes
[coeff,~,~] = pca(spoint);

meanXYZ = mean(spoint);
temp = spoint-meanXYZ;
temp_dist = sqrt(temp(:,1).^2+temp(:,2).^2+temp(:,3).^2);
temp_list = temp_dist < median(temp_dist);
meanXYZ2 = mean(spoint(temp_list>0,:));
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

point5 = spoint*coeff;

%% Circle fitting
[~,temp_center] = CircleFit_NoiseRemoval(point5);
[center, coef, rad] = CircleFitting(point5,temp_center(:,1:2));

temp = point5(:,1:2)-repmat(center,size(point5,1),1);
temp_dist = sqrt(temp(:,1).^2+temp(:,2).^2);
mdist0 = mean(temp_dist);

if mdist0 > 1500
    [rad,coef,center] = CircleFitting_Alternative(point5);
end

point5_2 = o_point * coeff;
temp = point5_2(:,1:2)-repmat(center,size(point5_2,1),1);
rad3 = ComputeAngularCoord(temp);

%% Measure inclination of approximation planes
number = 3;
int_rad = min(rad):(max(rad)-min(rad))/number:max(rad);
sita2 = zeros(number,1);
rad2 = (min(rad)+(max(rad)-min(rad))/number/2:(max(rad)-min(rad))...
    /number:max(rad)-(max(rad)-min(rad))/number/2)';

for i = 1:number
    f1 = rad>=int_rad(1,i);
    f2 = rad< int_rad(1,i+1);
    f = f1.*f2;
    temp_coef = pca(point5(f>0,:));
    h = temp_coef(:,3);
    if h(3)<0
        h = -h;
    end

    temp_rad = rad2(i);
    temp_R = [cos(temp_rad),sin(temp_rad)];

    dxdt = coef(2)*cos(temp_rad)+(coef(2)*temp_rad + coef(1))*(-sin(temp_rad));
    dydt = coef(2)*sin(temp_rad)+(coef(2)*temp_rad + coef(1))*(cos(temp_rad));
    temp_R2 = [dxdt dydt]/norm([dxdt dydt]);

    temp_R3 = [-temp_R2(2) temp_R2(1)];
    if dot(temp_R,temp_R3) < 0
        temp_R3 = -temp_R3;
    end
    
    cos_alpha = dot(temp_R3,[h(1),h(2)])/norm([h(1),h(2)]);
    a = norm([h(1),h(2)])*cos_alpha;
    sita2(i,1) = atan(a/h(3));
end

%% Make two fitting curves

xres2 = 1;
zres2 = 1;
margin = 0.1;
strad = min(rad3)-(max(rad3)-min(rad3))*margin;
enrad = max(rad3)+(max(rad3)-min(rad3))*margin;
a = coef(2); b = coef(1);
fun4 = @(sita) sqrt(a^2 + (a*sita + b).^2);
fun5 = @(alpha) integral(fun4,strad,strad+alpha);
xlen = fun5(enrad-strad);

t_list = (strad:0.01:enrad)';
y_list = zeros(size(t_list,1),1);
for i = 1:size(t_list,1)
    y_list(i,1) = fun5(t_list(i,1)-strad);
end
p = polyfit(y_list,t_list,3);
fun_X = @(X) p(1)*(X.^3)+p(2)*(X.^2)+p(3)*X + p(4);

width = round(55/PIXEL_WIDTH);
depth = round(18/Z_STEP);

temp_im = zeros(round(xlen/xres2),width*2+1,depth*2+1);
temp_im(:,width+1,depth+1)=1;
xyz_list = GetCoordOfPositivePixels(temp_im);
centerline1 = ProjectLinearizedCoordToOriginal(xyz_list,width, depth, PIXEL_WIDTH, Z_STEP ...
    , xres2, zres2, coef, mean(point5(:,3)), rad2, sita2, center, fun_X, coeff);

temp_im = zeros(round(xlen/xres2),width*2+1,depth*2+1);
temp_im(:,1,depth+1)=1;
xyz_list = GetCoordOfPositivePixels(temp_im);
centerline2 = ProjectLinearizedCoordToOriginal(xyz_list,width, depth, PIXEL_WIDTH, Z_STEP ...
    , xres2, zres2, coef, mean(point5(:,3)), rad2, sita2, center, fun_X, coeff);

% Convert to Cartesian coordinates
centerline1 = centerline1 * diag([PIXEL_WIDTH,PIXEL_WIDTH,Z_STEP]);
centerline2 = centerline2 * diag([PIXEL_WIDTH,PIXEL_WIDTH,Z_STEP]);

% Create vector list indicating inclinations of each point
inclines = (centerline2-centerline1);
inclines = inclines./sqrt(sum(inclines.^2,2));

