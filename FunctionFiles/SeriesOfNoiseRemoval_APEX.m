function pointList = SeriesOfNoiseRemoval_APEX(pointList, circleCenter, counterClockWise)
%*****************************************************************************
%This function uses "pointList", Cartesian coordinates of points,
%"circleCenter", pole of polar coodinate system, and "counterClockWise",
%whether cochlear is clockwise or counterClockWise.
%*****************************************************************************

%% Noise removal by radial coordinates
[radialCoords,angularCoords] = ConvertToPolar(pointList, circleCenter);
mdist = (prctile(radialCoords,10)+prctile(radialCoords,90))/2;
dist2 = radialCoords - ones(size(radialCoords,1),1)*mdist;
f3 = (dist2 < 70).*(dist2 > -100);
pointList = pointList(f3>0,:);

%% Noise removal by single-linkage clustering on polar coordinates
gain = 3; interval = 18;
h_dist = (angularCoords-median(angularCoords)).*radialCoords;
v_dist = gain*(radialCoords-median(radialCoords));
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
    temp2 = radialCoords(T==i>0);
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

pointList = pointList((f1.*f2)>0,:);

%% Noise removal by single-linkage clustering on Cartesian coordinates
Y = pdist(pointList,'euclid');
Z = linkage(Y,'single');
T = cluster(Z,'cutoff',15,'criterion','distance');
clusterSize = zeros(max(T),1);
tmppoints = [];
for i =1:max(T)
    clusterSize(i,1) = sum(T==i);
end
f = clusterSize > 10;
idx = find(f);
for i =1:sum(f)
    tmppoints = [tmppoints; pointList(T==idx(i),:)];
end
pointList = tmppoints;

%% Noise removal by angular coordinates
[~,angularCoords] = ConvertToPolar(pointList, circleCenter);
pointCounts = size(angularCoords,1);
if counterClockWise == 1
    for i = 1:size(angularCoords,1)
        temp = angularCoords(i);
        f1 = angularCoords > (temp-0.5);
        f2 = angularCoords < (temp+0.05);
        pointCounts(i,1) = sum(f1.*f2);
    end
else
    for i = 1:size(angularCoords,1)
        temp = angularCoords(i);
        f1 = angularCoords > (temp-0.05);
        f2 = angularCoords < (temp+0.5);
        pointCounts(i,1) = sum(f1.*f2);
    end
end
pointList(pointCounts<20,:) = [];