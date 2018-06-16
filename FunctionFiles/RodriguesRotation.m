function xyz = RodriguesRotation(xyz,vector,sita)
%*****************************************************************************
%This function uses "xyz", xyz coordinates of points, "vector", rotation
%axis, and "sita" rotation angle. The output is "xyz", rotated coordinates
%of points.
%*****************************************************************************

%% Do Rodrigues' rotaion
xyz = xyz';
n = vector./norm(vector);

a = cos(sita) + (n(1)^2)*(1-cos(sita));
b = n(1)*n(2)*(1-cos(sita)) - n(3)*sin(sita);
c = n(1)*n(3)*(1-cos(sita)) + n(2)*sin(sita);
d = n(2)*n(1)*(1-cos(sita)) + n(3)*sin(sita);
e = cos(sita) + (n(2)^2)*(1-cos(sita));
f = n(2)*n(3)*(1-cos(sita)) - n(1)*sin(sita);
g = n(3)*n(1)*(1-cos(sita)) - n(2)*sin(sita);
h = n(3)*n(2)*(1-cos(sita)) + n(1)*sin(sita);
i = cos(sita) + (n(3)^2)*(1-cos(sita));

R = [a b c; d e f; g h i];
xyz = R * xyz;
xyz = xyz';