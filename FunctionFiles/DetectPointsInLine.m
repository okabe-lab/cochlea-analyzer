function [length,width,point2] = DetectPointsInLine(points,temp_center,gain,interval)
%*****************************************************************************
%This function uses "points", coordinates of points which may represent
%organ of corti, "temp_center", initial value of center coordinate for
%circle fitting. "gain", gain for radial coordinate, "interval", cutoff
%value of distance for clustering. The outputs are "length", selected 
%cluster length, "width", selected cluster width, "point2", points in
%selected cluster.
%*****************************************************************************

%% Clustering with minimum distance method on radial and angular coordinates
temp = points(:,1:2)-repmat(temp_center,size(points,1),1);
rad = ComputeAngularCoord(temp);
dist = sqrt(temp(:,1).^2+temp(:,2).^2);

nrad = (rad-mean(rad))/(max(rad)-min(rad));
ndist = gain*(dist-mean(dist))/(max(dist)-min(dist));
npoints = [nrad ndist];

Y = pdist(npoints,'euclid');
Z = linkage(Y,'single');
T = cluster(Z,'cutoff',interval,'criterion','distance');
candid_group2 = cell(max(T),1);
for i =1:max(T)
    candid_group2{i,1} = npoints(T==i>0,:);
end

%% Select long and narrow cluster
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

if sum(tleng > 0.85) > 0
    temp = tnum.*(tleng > 0.85);
    temp(temp==0,:)=inf;
    [~,idx] = min(temp);
else
    [~,idx] = max(tleng);
end
point2 = points(T==idx,:);

center = ComputeCircleCenter(point2, temp_center, 7);

temp = point2(:,1:2)-repmat(center,size(point2,1),1);
dist2 = sqrt(temp(:,1).^2+temp(:,2).^2);
rad2 = ComputeAngularCoord(temp);

rad33 = prctile(rad2,33);
rad67 = prctile(rad2,67);
f1 = rad2 <= rad33;
f2 = (rad2 > rad33).*(rad2 <= rad67);
f3 = rad2 > rad67;
dist2_1 = max(dist2(f1>0,:))-min(dist2(f1>0,:));
dist2_2 = max(dist2(f2>0,:))-min(dist2(f2>0,:));
dist2_3 = max(dist2(f3>0,:))-min(dist2(f3>0,:));

length = tleng(idx);
width = max(dist2)-min(dist2)+std([dist2_1 dist2_2 dist2_3]);
