function [point5,center,point6,width] = CircleFit_NoiseRemoval(point4)
%*****************************************************************************
%This function uses "point4", xyz coordinates of points. The outputs are 
%"point5", new coordinates of points, "center", center of fitted circle, 
%"point6", new coordinates of points with restricted criteria, "width", 
%spread of points in radial coorodinate.
%*****************************************************************************

%% Circle fitting
r = (360:20:800)';
rsize = size(r,1);
wid = zeros(rsize,1);
tcent = zeros(rsize,2);
for i = 1:rsize
    [wid(i,1),tcent(i,:)] = SimpleCircleFitEvaluation(point4,r(i,1));
end
[~,idx] = min(wid);
temp_center = tcent(idx,:);
[center,~,~,width] = CircleFitting(point4, temp_center);

%% Remove noise
temp = point4(:,1:2)-repmat(center,size(point4,1),1);
dist = sqrt(temp(:,1).^2+temp(:,2).^2);
mdist = (prctile(dist,10)+prctile(dist,90))/2;

dist2 = dist - ones(size(dist,1),1)*mdist;
f3 = (dist2 < 70).*(dist2 > -100);

point5 = point4(f3>0,:);
point5 =  NoiseRemovalInPolarCoord(point5,center,3,18);

point6 = [];
Y = pdist(point5,'euclid');
Z = linkage(Y,'single');
T = cluster(Z,'cutoff',15,'criterion','distance');
mem_num = zeros(max(T),1);
for i =1:max(T)
    mem_num(i,1) = sum(T==i);
end
f = mem_num > 10;
idx = find(f);
for i =1:sum(f)
    point6 = [point6; point5(T==idx(i)>0,:)];
end