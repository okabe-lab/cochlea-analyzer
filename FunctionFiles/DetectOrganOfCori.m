function [centerLine,inclines,points,intensity] = DetectOrganOfCori(originalIm,PIXEL_WIDTH,Z_STEP,CounterClockWise)
%*****************************************************************************
%This function uses "originalIm", image stack, and "PIXEL_WIDTH" and "Z_STEP", image
%resolution. The outputs are "compline", fitting curve along organ of corti,
%"vector", vectors indicating inclination of approximation planes, "points"
%, intensity peaks, and "intensity", intensity of each peak point.
%*****************************************************************************

%% Detect organ of corti by intensity peak detection
% Obtain regional maxes of pixel values in image stack
[pointList,backGroundLevel] = ObtainPointGroup(originalIm,PIXEL_WIDTH,Z_STEP);

% Single-linkage clustering
y = pdist(pointList,'euclid');
z = linkage(y,'single');
assignedClusterNos = cluster(z,'cutoff',25,'criterion','distance');

% Filter by clulster size
clusterSizes = zeros(max(assignedClusterNos),1);
for i =1:max(assignedClusterNos)
    clusterSizes(i,1) = sum(assignedClusterNos==i);
end
idx = find(clusterSizes > 200);
candidateClusters = cell(size(idx));
for i =1:size(candidateClusters,1)
    candidateClusters{i,1} = pointList(assignedClusterNos==idx(i),:);
end

% Select cluster correspond to organ of corti
if size(candidateClusters,1) > 1
    idx = FindClusterFromCorti(candidateClusters,originalIm,backGroundLevel,PIXEL_WIDTH,Z_STEP);
    pointList = candidateClusters{idx(1),1};
else
    pointList = candidateClusters{1,1};
end

if 0
    % Search undetected hair cells around extracted points
    pointList = ReObtainPointGroup(pointList,originalIm,PIXEL_WIDTH,Z_STEP);
    
    %% Remove noise
    point3 = RemoveSmallClusters(pointList, 24, 21);
    
    for ii = 1:5
        coef = pca(point3);
        point4 = point3*coef;
        zcenter = median(point4(:,3));
        squareD = squareform(pdist(point4(:,1:2))); %2ŽŸŒ³
        tpoint3 = point3;
        for i = 1:size(point4,1)
            temp = point4(i,:);
            temp2 = squareD(:,i);
            temp2(i,1) = inf;
            [minimum,idx] = min(temp2);
            if minimum < 3
                z1 = abs(temp(1,3)-zcenter);
                temp3 = point4(idx,:);
                z2 = abs(temp3(1,3)-zcenter);
                [~,idx2] = min([z1,z2]);
                if idx2 == 2
                    tpoint3(i,:) = point3(idx,:);
                end
            end
        end
        point3 = unique(tpoint3,'rows');
    end
    
    squareD = squareform(pdist(point3));
    tpoint3 = point3;
    for i = 1:size(point3,1)
        temp = round(point3(i,:)*diag([1/PIXEL_WIDTH, 1/PIXEL_WIDTH, 1/Z_STEP]));
        temp2 = squareD(:,i);
        temp2(temp2==0) = inf;
        [minimum,idx] = min(temp2);
        if minimum < 5
            intens1 = originalIm(temp(1),temp(2),temp(3));
            temp3 = round(point3(idx,:)*diag([1/PIXEL_WIDTH, 1/PIXEL_WIDTH, 1/Z_STEP]));
            intens2 = originalIm(temp3(1),temp3(2),temp3(3));
            [~,idx2] = max([intens1,intens2]);
            if idx2 == 2
                tpoint3(i,:) = temp3;
            end
        end
    end
    point3 = unique(tpoint3,'rows');
    ppoint2 = point3*diag([1/PIXEL_WIDTH, 1/PIXEL_WIDTH, 1/Z_STEP]);
    tempim = originalIm(:,:,round(min(ppoint2(:,3))):round(max(ppoint2(:,3))));
    maxim = max(tempim,[],3);
    maxim2 = medfilt2(maxim);
    
    %% Template matching by image of hair cell
    template = [166    166    175    175    175    143    150;...
        315    492    592    592    555    378    175;...
        525    687    945    945    776    718    412;...
        625    800    986   1005    986    800    669;...
        625    800    983   1012    966    800    669;...
        474    687    703    945    800    706    412;...
        250    380    528    582    582    328    234];
    w = floor(size(template,1)/2);
    corrCoefMat = normxcorr2(template,maxim2);
    corrCoefMat = corrCoefMat(w+1:end-w,w+1:end-w);
    
    noise_temp = [       60    61    59    52    58    55    53;...
        48    55    52    59    57    55    54;...
        64    62    59    69    53    56    58;...
        74    62    74   405   113    62    63;...
        69    56    77   103    73    65    75;...
        80    70    64    74    96    65    58;...
        68   100    85    82   136    95    65];
    
    C2 = normxcorr2(noise_temp,maxim);
    SE = strel('square',3);
    C3 = imdilate(C2,SE);
    
    corr = zeros(size(point3,2));
    for i = 1:size(point3,1)
        corr(i,1) = corrCoefMat(round(point3(i,1)),round(point3(i,2)));
        corr(i,2) = C3(round(point3(i,1)),round(point3(i,2)));
    end
    
    f = (corr(:,1)>0.1) + (corr(:,2)<0.6)>1;
    point4 = point3(f>0,:);
    
    point4 = RemoveSmallClusters(point4, 20, 21);
    
    point3_2 = point3;
    point3 = point4;
    
    pointList = point3*diag([PIXEL_WIDTH, PIXEL_WIDTH, Z_STEP]);
    pointList2 = point3_2*diag([PIXEL_WIDTH, PIXEL_WIDTH, Z_STEP]);
    
    %% PCA and Circle fitting (determine center of circle)
    [pointList_PCA,pcCoefficients] = ProjectOntoPCA_Twice(pointList);
    
    zcenter = median(pointList_PCA(:,3));
    fz1 = pointList_PCA(:,3) > zcenter+20;
    fz2 = pointList_PCA(:,3) < zcenter-20;
    fz = fz1 + fz2;
    pointList_PCA(fz>0,:) = [];
    
    [pointList_PCA,circleCenter] = CircleFit_NoiseRemoval(pointList_PCA);
    
    %% Remove Noise
    point5_2 = pointList2;
    corr = zeros(size(point5_2,1),1);
    ppoint5_2 = point5_2*diag([1/PIXEL_WIDTH 1/PIXEL_WIDTH 1/Z_STEP]);
    for i = 1:size(point5_2,1)
        corr(i,1) = corrCoefMat(round(ppoint5_2(i,1)),round(ppoint5_2(i,2)));
    end
    
    point5_3 = pointList2*pcCoefficients;
    temp = point5_3(:,1:2)-repmat(circleCenter,size(point5_3,1),1);
    dist = sqrt(temp(:,1).^2+temp(:,2).^2);
    
    temp2 = pointList_PCA(:,1:2)-repmat(circleCenter,size(pointList_PCA,1),1);
    dist_2 = sqrt(temp2(:,1).^2+temp2(:,2).^2);
    
    distnum = floor((max(dist_2)-min(dist_2))/10);
    startdist = max(dist_2)-distnum*10;
    sumcor = zeros(distnum,1);
    numlist = zeros(distnum,1);
    distlist = zeros(distnum,1);
    for i = 1:distnum
        f1 = dist >= ((i-1)*10+startdist);
        f2 = dist <= (i*10+startdist);
        sumcor(i,1) = sum(corr(f1.*f2>0,:));
        numlist(i,1) = sum(f1.*f2);
        distlist(i,1) = i*10+startdist;
    end
    
    rate = sumcor./numlist;
    f1 = numlist > mean(numlist)*0.8;
    rate2 = rate.*f1;
    [~,idx] = sort(rate2,'descend');
    distth = max(distlist(idx(1:2),:))+25;
    pointList_PCA(dist_2 > distth,:) = [];
    
    %% Circle fitting (determine radius of circle)
    [dist,rad] = ConvertToPolar(pointList_PCA,circleCenter);
    
    number = 6;
    int_rad = min(rad):(max(rad)-min(rad))/number:max(rad);
    m_data = zeros(number,1);
    for i = 1:size(int_rad,2)-1
        f1 = rad>=int_rad(1,i);
        f2 = rad< int_rad(1,i+1);
        f = f1.*f2;
        temp = dist(f>0);
        m_data(i,1) = (prctile(temp,90) + prctile(temp,10))/2;
    end
    coef_g0 = [mean(m_data) 0]';
    
    mdist0 = mean(dist);
    if mdist0 > 1500
        [rad,coef_g0,circleCenter,pointList_PCA] = CircleFitting_Alternative(pointList_PCA);
    end
    
    %% Get intensities of peak points
    Isize = size(originalIm);
    ppoint5 = (pointList_PCA/pcCoefficients)*diag([1/PIXEL_WIDTH 1/PIXEL_WIDTH 1/Z_STEP]);
    temppoint = round(ppoint5);
    temppoint = sub2ind(Isize,temppoint(:,1),temppoint(:,2),temppoint(:,3));
    intlist = originalIm(temppoint);
    intensity = prctile(intlist,90);
    
end

Isize = size(originalIm);
im0 = originalIm;
[point3, C, point3_2] = ReObtainPointGroup2(pointList,originalIm,PIXEL_WIDTH,Z_STEP,0);

[point4,coeff] = urata_newf13(point3);
zcenter = median(point4(:,3));
fz1 = point4(:,3) > zcenter+20;
fz2 = point4(:,3) < zcenter-20;
fz = fz1 + fz2;
point4(fz>0,:) = [];

[point5,center_g0,~,width] = urata_newf68(point4,0);
temp = point5(:,1:2)-repmat(center_g0,size(point5,1),1);
rad = calc_rad(temp);
dist = sqrt(temp(:,1).^2+temp(:,2).^2);

if width > 100 || max(rad)-min(rad) < 0.5 || median(dist) < 300 || median(dist) > 1500
    [width median(dist)]
    point4_2 = [-point4(:,1) -point4(:,2) point4(:,3)];
    
    r = (360:20:800)';
    rsize = size(r,1);
    wid = zeros(rsize,1);
    tcent = zeros(rsize,2);
    for i = 1:rsize
        [wid(i,1),tcent(i,:)] = urata_newf81(point4,r(i,1));
    end
    
    [~,idx] = min(wid);
    temp_center = tcent(idx,:);
    
    point4_3 =  urata_newf76(point4_2,temp_center,3,20,0);
    
    [point5_2,center_g0_2,~,width_2] = urata_newf68(point4_3,0);
    temp = point5_2(:,1:2)-repmat(center_g0_2,size(point5_2,1),1);
    dist_2 = sqrt(temp(:,1).^2+temp(:,2).^2);
    if width_2 < 100 && median(dist_2) > 300 && median(dist_2) < 1500
        point5 = point5_2;
        center_g0 = center_g0_2;
        coeff = [-coeff(:,1) -coeff(:,2) coeff(:,3)];
    end
end

point2_2 = pointList*coeff;
temp = point2_2(:,1:2)-repmat(center_g0,size(point2_2,1),1);
dist1 = sqrt(temp(:,1).^2+temp(:,2).^2);

temp2 = point5(:,1:2)-repmat(center_g0,size(point5,1),1);
rad2 = calc_rad(temp2);
dist2 = sqrt(temp2(:,1).^2+temp2(:,2).^2);

f1 = dist2 >= min(dist1);
if max(dist1)-min(dist1) > 80
    f2 = dist2 <= max(dist1);
else
    f2 = dist2 <= min(dist1)+100;
end
point5(f1.*f2==0,:) = [];

point5_2 = point3_2;
corr = zeros(size(point5_2,1),1);
ppoint5_2 = point5_2*diag([1/PIXEL_WIDTH 1/PIXEL_WIDTH 1/Z_STEP]);
for i = 1:size(point5_2,1)
    corr(i,1) = C(round(ppoint5_2(i,1)),round(ppoint5_2(i,2)));
end
f = corr>0.6;

point5_3 = point3_2*coeff;
temp = point5_3(:,1:2)-repmat(center_g0,size(point5_3,1),1);
rad = calc_rad(temp);
dist = sqrt(temp(:,1).^2+temp(:,2).^2);

temp2 = point5(:,1:2)-repmat(center_g0,size(point5,1),1);
rad_2 = calc_rad(temp2);
dist_2 = sqrt(temp2(:,1).^2+temp2(:,2).^2);

distnum = floor((max(dist_2)-min(dist_2))/10);
startdist = max(dist_2)-distnum*10;
sumcor = zeros(distnum,1);
numlist = zeros(distnum,1);
distlist = zeros(distnum,1);
for i = 1:distnum
    f1 = dist >= ((i-1)*10+startdist);
    f2 = dist <= (i*10+startdist);
    sumcor(i,1) = sum(corr(f1.*f2>0,:));
    numlist(i,1) = sum(f1.*f2);
    distlist(i,1) = i*10+startdist;
end

rate = sumcor./numlist;
f1 = numlist > mean(numlist)*0.8;
f2 = sumcor./numlist > mean(rate(f1>0,:));
rate2 = rate.*f1;
[~,idx] = sort(rate2,'descend');
f3 = zeros(size(f1));
f3(idx(1:2),:) = 1;
distth = max(distlist(idx(1:2),:))+25;

point5(dist_2 > distth,:) = [];

%% Circle Fitting
temp = point5(:,1:2)-repmat(center_g0,size(point5,1),1);
rad = calc_rad(temp);
dist = sqrt(temp(:,1).^2+temp(:,2).^2);
number = 6;
int_rad = min(rad):(max(rad)-min(rad))/number:max(rad);
m_data = zeros(number,1);
for i = 1:size(int_rad,2)-1
    f1 = rad>=int_rad(1,i);
    f2 = rad< int_rad(1,i+1);
    f = f1.*f2;
    temp = dist(f>0);
    m_data(i,1) = (prctile(temp,90) + prctile(temp,10))/2;
end

coef_g0 = [mean(m_data) 0]';

temp = point5(:,1:2)-repmat(center_g0,size(point5,1),1);
temp_dist = sqrt(temp(:,1).^2+temp(:,2).^2);
mdist0 = mean(temp_dist);

if mdist0 > 1500
    [rad,coef_g0,center_g0,point5] = urata_newf59(point4);
end

%% Get brightness
ppoint5 = (point5/coeff)*diag([1/PIXEL_WIDTH 1/PIXEL_WIDTH 1/Z_STEP]);
temppoint = round(ppoint5);
temppoint = sub2ind(Isize,temppoint(:,1),temppoint(:,2),temppoint(:,3));
intlist = im0(temppoint);
intensity = prctile(intlist,90);

pointList_PCA = point5;
circleCenter = center_g0;
pcCoefficients = coeff;

%% Remove Noise
f_g0 = IsNearBoundaries(pointList_PCA,circleCenter,coef_g0,pcCoefficients,Isize*diag([PIXEL_WIDTH,PIXEL_WIDTH,Z_STEP]));
temppoint = pointList_PCA(f_g0==0,:);
Y = pdist(temppoint,'euclid');
Z = linkage(Y,'single');
T = cluster(Z,'cutoff',60,'criterion','distance');
if max(T)>1
    mem_num = zeros(max(T),1);
    for i =1:max(T)
        mem_num(i,1) = sum(T==i);
    end
    [~,idx] = max(mem_num);
    temppoint = temppoint(T==idx,:);
    [~,temp] = knnsearch(temppoint, pointList_PCA);
    f_g0 = temp~=0;
end

%% Divide points into several groups
radw = max(rad(f_g0==0,:))-min(rad(f_g0==0,:));
radth1 = 0.6;
radth2 = 1.5;
if radw < radth1
    p1 = 0;
    group1 = pointList_PCA;
    group2 = pointList_PCA;
    group1_2 = pointList_PCA((f_g0 == 0)>0,:);
    group2_2 = pointList_PCA((f_g0 == 0)>0,:);
else
    if radw > radth2
        p1 = 42;
    else
        p1 = 42*(radw-radth1)/(radth2-radth1);
    end
    p2 = 100-p1;
    rad40 = prctile(rad(f_g0==0,:),p1);
    rad60 = prctile(rad(f_g0==0,:),p2);
    group1 = pointList_PCA(rad < rad60,:);
    group2 = pointList_PCA(rad > rad40,:);
    group1_2 = pointList_PCA((rad < rad60).*(f_g0 == 0)>0,:);
    group2_2 = pointList_PCA((rad > rad40).*(f_g0 == 0)>0,:);
end

if p1 == 0
    [group1_3,center_g1] = CircleFit_NoiseRemoval2(group1_2,circleCenter);
    [center_g1, coef_g1] = CircleFitting(group1_3,center_g1);
    f_g1 = IsNearBoundaries(group1_3,center_g1,coef_g1,pcCoefficients,Isize*diag([PIXEL_WIDTH,PIXEL_WIDTH,Z_STEP]));
    group1_4 = group1_3(f_g1==0,:);
    group2_4 = group1_4;
else
    [group1_3,center_g1] = CircleFit_NoiseRemoval2(group1_2,circleCenter);
    [center_g1, coef_g1] = CircleFitting(group1_3,center_g1);
    f_g1 = IsNearBoundaries(group1_3,center_g1,coef_g1,pcCoefficients,Isize*diag([PIXEL_WIDTH,PIXEL_WIDTH,Z_STEP]));
    group1_4 = group1_3(f_g1==0,:);
    
    [group2_3,center_g2] = CircleFit_NoiseRemoval2(group2_2,circleCenter);
    [center_g2, coef_g2] = CircleFitting(group2_3,center_g2);
    f_g2 = IsNearBoundaries(group2_3,center_g2,coef_g2,pcCoefficients,Isize*diag([PIXEL_WIDTH,PIXEL_WIDTH,Z_STEP]));
    group2_4 = group2_3(f_g2==0,:);
end

group1 = group1 / pcCoefficients;
selectedGroup1 = group1_4 / pcCoefficients;
group2 = group2 / pcCoefficients;
selectedGroup2 = group2_4 / pcCoefficients;

%% Make fitting curves for each group of peaks
if p1 == 0
    [centerLine,inclines] = DrawLinesAlongCorti(group1,selectedGroup1,PIXEL_WIDTH,Z_STEP);
else
    [centerLine1,inclines1] = DrawLinesAlongCorti(group1,selectedGroup1,PIXEL_WIDTH,Z_STEP);
    [centerLine2,inclines2] = DrawLinesAlongCorti(group2,selectedGroup2,PIXEL_WIDTH,Z_STEP);
    
    [centerLine,inclines] = MergeLines(centerLine1,centerLine2,inclines1,inclines2);
end

if CounterClockWise == 0
    centerLine = flipud(centerLine);
    inclines = flipud(inclines);
end

points = pointList_PCA / pcCoefficients;
