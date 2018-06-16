function [rss,center_x,center_y,angular_coord,radial_coord] ...
    = ComputeRssOfCochlearHelixFit(point_group,rotation_rad1,rotation_rad2)
%*****************************************************************************
%This function uses "point_group", cartesian coordinate of helically
%distributed points, "rotation_rad1", first rotation angle, "rotation_rad2"
%second rotation angle. The outputs are "rss" residual sum of squares of
%fitting, and "angular_coord" and "radial_coord", polar coordinates of point
%group. 
%*****************************************************************************

%% Rotate point_group by specified angles
temp = RodriguesRotation(point_group,[1,0,0],rotation_rad1);
rot_xyz = RodriguesRotation(temp,[0,1,0],rotation_rad2);

%% Circle fitting
X3 = rot_xyz(:,1:2);
x = X3(:,1);
y = X3(:,2);
a0 = [mean(x),mean(y),max(x)-mean(x)];
f = @(a) norm((x-a(1)).^2 + (y-a(2)).^2 - a(3).^2); 
af1 = fminsearch(f, a0); 

rot_center = [af1(1),af1(2)];
center_x = af1(1);
center_y = af1(2);

xyz3 = [1 0 -rot_center(1,1); 0 1 -rot_center(1,2); 0 0 1]*[X3 ones(size(X3,1),1)]';
xyz3 = xyz3(1:2,:)';

%% Transform cartesian coodinates into polar coordinates

angular_coord = zeros(size(X3,1),1);
sita2 = zeros(size(X3,1),1);
radial_coord = zeros(size(X3,1),1);
for i = 1:size(X3,1)
    temp = xyz3(i,1:2);
    radial_coord(i,1) = norm(temp);
    
    temp = temp/norm(temp);
    
    if i == 1
        v0 = temp;
        v00 = temp;
        angular_coord(i,1) = 0;
        sita2(i,1) = 0;
        judge = cross([v00 0],[temp 0]);
        continue
    end

    temp2 = acos(dot(v0,temp));
    v0 = temp;
    angular_coord(i,1) = angular_coord(i-1,1) + temp2;

    temp3 = acos(dot(v00,temp));
    temp4 = cross([v00 0],[temp 0]);
    
    if judge(3) > 0
        if temp4(1,3)<0
            temp3 = 2*pi-temp3;
        end
    else
        if temp4(1,3)>0
            temp3 = 2*pi-temp3;
        end
    end
    sita2(i,1) = temp3;
end

%% Calculate RSS of fitting
ten = round(size(angular_coord,1)/50);
sita0 = angular_coord(ten:end-ten,:);
dist0 = radial_coord(ten:end-ten,:);

X = [ones(length(sita0),1) sita0];
b = X\dist0;
yCalc = X*b;

yresid = dist0-yCalc;
rss = sum(yresid.^2);