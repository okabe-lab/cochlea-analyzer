function newPointList = RemoveSmallClusters(pointList, cutOffDist, cutOffSize)
%*****************************************************************************
%This function uses "pointList", Cartesian coordinates of points,
%"cutOffDist", cut off distance for single-linkage clustering, and "cutOffSize"
%, minimum cluster size to be kept.  
%*****************************************************************************

% Single-linkage clustering
y = pdist(pointList,'euclid');
z = linkage(y,'single');
assignedClusterNos = cluster(z,'cutoff',cutOffDist,'criterion','distance');

clusterSizes = zeros(max(assignedClusterNos),1);
for i =1:max(assignedClusterNos)
    clusterSizes(i,1) = sum(assignedClusterNos==i); % Check cluster sizes
end
idx = find(clusterSizes >= cutOffSize); % Find clusters to be kept
newPointList = [];
for i =1:numel(idx)
    newPointList = [newPointList; pointList(assignedClusterNos==idx(i),:)];
end