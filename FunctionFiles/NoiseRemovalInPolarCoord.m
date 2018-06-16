function point2 = NoiseRemovalInPolarCoord(points,temp_center,gain,interval)
%*****************************************************************************
%This function uses "points", coordinates of points which may represent 
%organ of corti, "temp_center", coordinates of center of fitted circle,
%"gain", gain for radial coordinate, and "interval", cutoff value of 
%distance for clustering. The output is "point2", coordinates of selected
%points.
%*****************************************************************************

%Remove noise by clustering with minimum distance method
%on radial and angular coordinates
temp = points(:,1:2)-repmat(temp_center,size(points,1),1);
rad = ComputeAngularCoord(temp);
dist = sqrt(temp(:,1).^2+temp(:,2).^2);
h_dist = (rad-median(rad)).*dist;
v_dist = gain*(dist-median(dist));
npoints = [h_dist v_dist];

Y = pdist(npoints,'euclid');
Z = linkage(Y,'single');
T = cluster(Z,'cutoff',interval,'criterion','distance');
candid_group2 = cell(max(T),1);
for i =1:max(T)
    candid_group2{i,1} = npoints(T==i>0,:);
end

tleng = zeros(size(candid_group2,1),1);
tnum = zeros(size(candid_group2,1),1);
thick = zeros(size(candid_group2,1),1);
for i = 1:size(candid_group2,1)
    temp = candid_group2{i,1};
    tleng(i,1) = max(temp(:,1))-min(temp(:,1));
    tnum(i,1) = size(temp,1);
    temp2 = dist(T==i>0);
    thick(i,1) = max(temp2)-min(temp2);
end

point2_2 = [];
temp = tleng > max(tleng)*0.3;
if sum(temp) > 0
    for i = 1:size(candid_group2,1)
        if temp(i,1) == 1
            point2_2 = [point2_2; candid_group2{i,1}];
        end
    end
    f1 = v_dist <= max(point2_2(:,2));
    f2 = v_dist >= min(point2_2(:,2))-50;
else
    f1 = ones(size(v_dist));
    f2 = ones(size(v_dist));
end

point2 = points((f1.*f2)>0,:);



