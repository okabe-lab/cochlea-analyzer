function [centerLine,inclines,pointList,imBrightness] ...
    = DetectOrganOfCori_Apex(originalIm,PIXEL_WIDTH,Z_STEP,counterClockWise)
%*****************************************************************************
%This function uses "originalIm", image stack, "PIXEL_WIDTH" and "Z_STEP", 
%image resolution, and "counterClockWise", whether cochlear is clockwise or
%counterClockWise. The outputs are "centerLine", fitting curve along organ of
%corti, "inclines", vectors indicating inclinations of approximation planes, 
%"pointList" extracted points, and "imBrightness", brightness of image stack.
%This fuction is made for the image stack taken from apical end.
%*****************************************************************************

%% Detect organ of corti
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

%% Remove noise with template matching
% Create maximum intensity projection image and resize
imSize = size(originalIm);
zBegin = round(min(pointList(:,3))/Z_STEP);
zEnd = round(max(pointList(:,3))/Z_STEP);
mipIm = medfilt2(max(originalIm(:,:,zBegin:zEnd),[],3)); % MIP + Median filtering
mipIm = imresize(mipIm,ceil([imSize(1)/PIXEL_WIDTH imSize(2)/PIXEL_WIDTH]));

% Template matching by image of typical outer hair cell
templateIm = [
    166    166    175    175    175    143    150;...
    315    492    592    592    555    378    175;...
    525    687    945    945    776    718    412;...
    625    800    986   1005    986    800    669;...
    625    800    983   1012    966    800    669;...
    474    687    703    945    800    706    412;...
    250    380    528    582    582    328    234];
coefficientsMat = TemplateMatching2D(mipIm, templateIm);

% Convert subscripts to linear index for collective acquisition of correlation coefficents
idx = sub2ind([imSize(1),imSize(2)],round(pointList(:,1)),round(pointList(:,2)));
coefficients = coefficientsMat(idx);
pointList = pointList(coefficients > 0,:); % Remove points with low correltion coefficients

% Search close pair of points and remove one of them
distList = squareform(pdist(pointList(:,1:2))); % Compute pairwise distances in 2D
pixList = round(pointList*diag([1/PIXEL_WIDTH 1/PIXEL_WIDTH 1/Z_STEP])); % Convert to pixel positions
for i = 1:size(pixList,1)
    dists = distList(:,i);
    dists(i) = inf; % Omit oneself
    [minDist,idx] = min(dists); % Search closest point
    if minDist < 5 % Select close pair within threshold of distance
        intens1 = originalIm(pixList(i,1),pixList(i,2),pixList(i,3));
        intens2 = originalIm(pixList(idx,1),pixList(idx,2),pixList(idx,3));
        if intens2 > intens1 % Replace when the other point is brighter
            pointList(i,:) = pointList(idx,:);
        end
    end
end
pointList = unique(pointList,'rows'); % Remove duplicates

%% Project onto PCA space 
[pointList_PCA,pcCoefficients] = ProjectOntoPCA(pointList);

% Circle fitting 
initialCenter = ComputeIntialCenter(pointList_PCA);
circleCenter = CircleFitting(pointList_PCA, initialCenter);

% Remove outliers in polar coordinates
pointList_PCA = SeriesOfNoiseRemoval_APEX(pointList_PCA, circleCenter, counterClockWise);

% Convert to polar coordinates
[radialCoords,angularCoords] = ConvertToPolar(pointList_PCA, circleCenter);

% If circle fitting is failed 
if mean(radialCoords) > 1500 % Judge by radius of fitted circle (usually < 1000 micron)
    [angularCoords,~,circleCenter,pointList_PCA] = CircleFitting_Alternative(pointList_PCA);
end

%% Get pixel values of extracted points
pointList = pointList_PCA / pcCoefficients; % Back to original coordinates and pixel positions
pointList_pixel = round((pointList)*diag([1/PIXEL_WIDTH, 1/PIXEL_WIDTH, 1/Z_STEP])); 
idx = sub2ind(imSize,pointList_pixel(:,1),pointList_pixel(:,2),pointList_pixel(:,3));
intensityList = originalIm(idx); % Get pixel values
imBrightness = prctile(intensityList,90); % Consider 90th percentile of values as brightness of image

%% Divide points into several groups
% Decide number of groups by range of 1st principal component (from 1 to 3)
pc1Range = max(pointList_PCA(:,1))-min(pointList_PCA(:,1)); %Range on 1st principal component
if pc1Range < 200 % One group
    group1 = pointList_PCA / pcCoefficients; % Back to original Cartesian coordinates
    selectedGroup1 = SelectPointsForFit(pointList_PCA, circleCenter, pcCoefficients, imSize ...
        , PIXEL_WIDTH, Z_STEP); % Remove points near boundaries of image stack
    
    groupCount = 1;

elseif pc1Range < 480 % Divide points into two groups
    % Set rate of margin for division
    angleRange = max(angularCoords)-min(angularCoords); 
    if angleRange > 1.5
        marginRate = 42;
    else
        marginRate = 42*(angleRange-0.6)/0.9;
    end
    angleTh1 = prctile(angularCoords,marginRate);
    angleTh2 = prctile(angularCoords,100 - marginRate);
    
    g1_PCA = pointList_PCA(angularCoords < angleTh2,:); % Group1
    group1 = g1_PCA / pcCoefficients; % Back to original Cartesian coordinates
    selectedGroup1 = SelectPointsForFit(g1_PCA, circleCenter, pcCoefficients, imSize ...
        , PIXEL_WIDTH, Z_STEP); % Remove points near boundaries of image stack
    
    g2_PCA = pointList_PCA(angularCoords > angleTh1,:); % Group2
    group2 = g2_PCA / pcCoefficients;
    selectedGroup2 = SelectPointsForFit(g2_PCA, circleCenter, pcCoefficients, imSize ...
        , PIXEL_WIDTH, Z_STEP);
    
    groupCount = 2;
    
else % Divide points into three groups
    angleTh30 = prctile(angularCoords,30); angleTh40 = prctile(angularCoords,40);
    angleTh60 = prctile(angularCoords,60); angleTh70 = prctile(angularCoords,70);

    g1_PCA = pointList_PCA(angularCoords < angleTh40,:); % Group1
    group1 = g1_PCA / pcCoefficients; % Back to original Cartesian coordinates
    selectedGroup1 = SelectPointsForFit(g1_PCA, circleCenter, pcCoefficients, imSize ...
        , PIXEL_WIDTH, Z_STEP); % Remove points near boundaries of image stack

    g2_PCA = pointList_PCA((angularCoords > angleTh30).*(angularCoords < angleTh70)>0,:); % Group2
    group2 = g2_PCA / pcCoefficients;
    selectedGroup2 = SelectPointsForFit(g2_PCA, circleCenter, pcCoefficients, imSize ...
        , PIXEL_WIDTH, Z_STEP);
        
    g3_PCA = pointList_PCA(angularCoords > angleTh60,:); % Group3
    group3 = g3_PCA / pcCoefficients;
    selectedGroup3 = SelectPointsForFit(g3_PCA, circleCenter, pcCoefficients, imSize ...
        , PIXEL_WIDTH, Z_STEP);
    
    groupCount = 3;
end

%% Make fitting curves for each group of peaks
if groupCount == 1
    [centerLine,inclines] = DrawLinesAlongCorti(group1,selectedGroup1,PIXEL_WIDTH,Z_STEP);
    
elseif groupCount == 2
    [centerLine1,inclines1] = DrawLinesAlongCorti(group1,selectedGroup1,PIXEL_WIDTH,Z_STEP);
    [centerLine2,inclines2] = DrawLinesAlongCorti(group2,selectedGroup2,PIXEL_WIDTH,Z_STEP);
    
    [centerLine,inclines] = MergeLines_Apex(centerLine1,centerLine2,inclines1,inclines2);
    
else
    [centerLine1,inclines1] = DrawLinesAlongCorti(group1,selectedGroup1,PIXEL_WIDTH,Z_STEP);
    [centerLine2,inclines2] = DrawLinesAlongCorti(group2,selectedGroup2,PIXEL_WIDTH,Z_STEP);
    [centerLine3,inclines3] = DrawLinesAlongCorti(group3,selectedGroup3,PIXEL_WIDTH,Z_STEP);
    
    [centerLine12,inclines12] = MergeLines_Apex(centerLine1,centerLine2,inclines1,inclines2);
    [centerLine,inclines] = MergeLines_Apex(centerLine12,centerLine3,inclines12,inclines3);
end

if counterClockWise == 0 
    centerLine = flipud(centerLine);
    inclines = flipud(inclines);
end

%% Adjustment for rapid change of curvature around apical end
adjustRange = 120; adjustZ = 15; adjustProx = 17;

% Gradually reduce z value of center line around apical end
centerLine(1:adjustRange,3)= centerLine(1:adjustRange,3)...
    - interp1([1 adjustRange],[adjustZ 1],(1:adjustRange)'); 
centerLine(:,3) = movmean(centerLine(:,3),3);

% Gradually shift center line toward the proximal direction
proxDirect = inclines(1,1:2)/norm(inclines(1,1:2)); % Direction vector toward modiolus
adjustVector = adjustProx * proxDirect; % Direction vector with length of "adjustProx"
centerLine(1:adjustRange,1:2) = centerLine(1:adjustRange,1:2) ...
    + interp1([1,adjustRange],[adjustVector;[0 0]],(1:adjustRange)');

% Adjust spacing of points on center line
cumulativeDist = [0; cumsum(sqrt(sum(diff(centerLine).^2,2)))];
cumulativeDist = cumulativeDist*(size(centerLine,1)/cumulativeDist(end));
centerLine = interp1(cumulativeDist,centerLine,(1:size(centerLine,1))');

% Adjust inclination vectors
inclines = interp1(cumulativeDist,inclines,(1:size(centerLine,1))');
inclines = movmean(inclines,5);

