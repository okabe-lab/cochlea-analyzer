function rad = ComputeAngularCoord(vector)
%*****************************************************************************
%This function uses "vector", xy coordinates of points, and the output is
%"rad", angular coordinate of polar coordinate system.
%*****************************************************************************

%% Transform xy coordinates into polar coordinates
vector = vector./sqrt(vector(:,1).^2+vector(:,2).^2);
a = acos(vector(:,1));
f1 = vector(:,2)<0;
a(f1>0) = 2*pi - a(f1>0);

if size(a,1)>2
    if max(a) - min(a) > 1.5 * pi
        f2 = (a < pi/2).*(a >= 0);
        a(f2>0) = 2*pi + a(f2>0);
    end
end
rad = a;
