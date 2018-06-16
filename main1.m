%*****************************************************************************
%This script create 1st linearized image of organ of corti ("linearizedIm1.tif").
%The script uses image stacks taken from cleared cochlear, and an excel
%file containing imaging positions on XY scanning stage of each image stack.
%The script also create some MAT-Files, which will use in succeeding scripts.
%*****************************************************************************

if not(exist('mainPath','var'))
    mainPath = fileparts(mfilename('fullpath')); % Get fullpath of this script
end
cd(mainPath)
addpath('.\FunctionFiles') % Add folder of fuction files to search path

%% Information input
% Specify the full path of a ImFolderPath containing image stacks and an excel file
% with coordinates for each image stack on XY scanning stage.
ImFolderPath = '.\TestData';
ExcelFileName = 'xyz.xlsx';
PIXEL_WIDTH = 1; % image resolution (pixel width, micron)
Z_STEP = 1; % image resolution (z step, micron)

% Adjustment for mismatch of refraction index between water and tissue-clearing reagent
Z_STEP = Z_STEP * (1.47/1.33);

% Make new folder
SubFolderPath = [ImFolderPath '\Result'];
if ~exist(SubFolderPath,'file')
    mkdir(ImFolderPath, '\Result')
end

%% Determine wether the cochlear is clockwise or counterClockWise
[imFileNames, stageCoord] = ObtainImStacksInfo(ImFolderPath,ExcelFileName);
f = FindStringPattern(imFileNames,"linearizedIm1");
imFileNames(f>0,:) = [];

% Get shifts between image stacks
shifts = zeros(size(imFileNames,1)-1,2);
for i = 1:size(imFileNames,1)-1
    shifts(i,:) = stageCoord(i,:) - stageCoord(i+1,:);
end
% Convert to pixel position
shiftsAtXYStage = round(shifts(:,1:2)*diag([1/PIXEL_WIDTH 1/PIXEL_WIDTH]));

% Judge clockwise or counterClockWise by calculating cross product
temp1 = [shiftsAtXYStage(1,:) 0]; temp2 = [shiftsAtXYStage(2,:) 0];
temp3 = [shiftsAtXYStage(3,:) 0]; temp4 = [shiftsAtXYStage(4,:) 0];
temp = cross(temp1,temp2) + cross(temp2,temp3) + cross(temp3,temp4);
counterClockWise = temp(3) < 0;

%% Detect organ of corti in each image stack

% Obtain curved lines along hair cells, point cloud of hair cells and
% brightness of each image stack.
imSizeList = zeros(size(imFileNames,1),3);
imBrightnessList = zeros(size(imFileNames,1),1);
for i = 1:size(imFileNames,1)
    fprintf('Detecting organ of Corti in image stack [%d / %d] ...\n',i,numel(imFileNames))

    originalIm = ImRead3D([ImFolderPath '\' imFileNames{i}]);
    imSizeList(i,:) = size(originalIm);
    
    if i == 1 % For image stack taken from apical end
        [centerLinePerIm,inclinesPerIm,pointList,imBrightnessList(i)] ...
            = DetectOrganOfCori_Apex(originalIm,PIXEL_WIDTH,Z_STEP,counterClockWise);
    else
        [centerLinePerIm,inclinesPerIm,pointList,imBrightnessList(i)] ...
            = DetectOrganOfCori(originalIm,PIXEL_WIDTH,Z_STEP,counterClockWise);
    end
    save([SubFolderPath '\' TrimTif(imFileNames{i}) '.mat'],'originalIm','pointList' ...
        ,'centerLinePerIm','inclinesPerIm','imSizeList','imBrightnessList','-v7.3');
end

%% Calculate accurate image shifts
if not(exist('shiftOfIms_pixel','var'))
    % Record accumulated shift values from first image stacks
    shiftOfIms_pixel = zeros(size(imFileNames,1),3);
end

% Load previous image stack, extracted points corresponding to hair cells.
load([SubFolderPath '\' TrimTif(imFileNames{1,1}) '.mat'],'-mat','originalIm','pointList');
prevIm = originalIm;
prevPoints = pointList * diag([1/PIXEL_WIDTH 1/PIXEL_WIDTH 1/Z_STEP]); % Convert to pixel position
prevImSize = size(prevIm);

% Compute shift between previous and current image stacks by cross-correlation.
for i = 2:size(imFileNames,1)
    fprintf('Calculating shifts between image stacks [%d and %d] ...\n',i-1,i)
    
    % Load image stack, extracted points corresponding to hair cells.
    load([SubFolderPath '\' TrimTif(imFileNames{i,1}) '.mat'],'-mat','originalIm','pointList');
    currentIm = originalIm;
    currentPoints = pointList * diag([1/PIXEL_WIDTH 1/PIXEL_WIDTH 1/Z_STEP]); % Convert to pixel position
    currentImSize = size(currentIm);
    
    % Estimate rough shift value of z coordinate from points extracted from image stacks
    % Select points from previous image whitin overlapped region
    shiftedPoints = prevPoints ...
        - repmat([shiftsAtXYStage(i-1,2) shiftsAtXYStage(i-1,1) 0],size(prevPoints,1),1);
    isInOverlap1 = (shiftedPoints(:,1)>0) .* (shiftedPoints(:,1)<currentImSize(1));
    isInOverlap2 = (shiftedPoints(:,2)>0) .* (shiftedPoints(:,2)<currentImSize(2));
    slectedPrevPoints = prevPoints(isInOverlap1.*isInOverlap2>0,:);
    shiftedPrevPoints = shiftedPoints(isInOverlap1.*isInOverlap2>0,:);
    
    % Select points from current image within overlapped region
    shiftedPoints = currentPoints ...
        + repmat([shiftsAtXYStage(i-1,2) shiftsAtXYStage(i-1,1) 0],size(currentPoints,1),1);
    isInOverlap1 = (shiftedPoints(:,1)>0) .* (shiftedPoints(:,1)<prevImSize(1));
    isInOverlap2 = (shiftedPoints(:,2)>0) .* (shiftedPoints(:,2)<prevImSize(2));
    selectedCurrentPoints = currentPoints(isInOverlap1.*isInOverlap2>0,:);
    
    % Removing outliers:
    % Compute distances between each point of previous image and the
    % corresponding closest point of current image.
    [~,minDist] = knnsearch(selectedCurrentPoints(:,1:2),shiftedPrevPoints(:,1:2));
    isInRange = minDist < 6; % Cut-off distance
    slectedPrevPoints(isInRange==0,:) = []; shiftedPrevPoints(isInRange==0,:) = [];
    % Vice versa
    [~,minDist] = knnsearch(shiftedPrevPoints(:,1:2),selectedCurrentPoints(:,1:2));
    isInRange = minDist < 6; % Cut-off distance
    selectedCurrentPoints(isInRange==0,:) = [];
    
    % Cut out image from overlapped region between two image stacks.
    overlapWidth1 = max(selectedCurrentPoints(:,1))-min(selectedCurrentPoints(:,1));
    overlapWidth2 = max(selectedCurrentPoints(:,2))-min(selectedCurrentPoints(:,2));
    if overlapWidth1 > 20 && overlapWidth2 > 20
        cutOutBegin1 = round(min(selectedCurrentPoints(:,1))+5);
        cutOutEnd1 = round(max(selectedCurrentPoints(:,1))-5);
        cutOutBegin2 = round(min(selectedCurrentPoints(:,2))+5);
        cutOutEnd2 = round(max(selectedCurrentPoints(:,2))-5);
    else
        % When little overlap is estimated between image stacks, use XY
        % scanning stage coordinates ("shiftsAtXYStage").
        cutOutBegin1 = round(max(1,1-shiftsAtXYStage(i-1,2)));
        cutOutEnd1 = round(min(prevImSize(1),prevImSize(1)-shiftsAtXYStage(i-1,2)));
        cutOutBegin2 = round(max(1,1-shiftsAtXYStage(i-1,1)));
        cutOutEnd2 = round(min(prevImSize(2),prevImSize(2)-shiftsAtXYStage(i-1,1)));
    end
    cutOut3 = round((min(selectedCurrentPoints(:,3))+max(selectedCurrentPoints(:,3)))/2);
    % Cut out 2D image from current image
    overlapIm = currentIm(cutOutBegin1:cutOutEnd1,cutOutBegin2:cutOutEnd2,cutOut3);
    
    % Compute shift of z coordinate between image stacks (provisional value)
    tempZShift = median(slectedPrevPoints(:,3))-median(selectedCurrentPoints(:,3));
    % Expected Z coordinate of "overlapIm" in previous image
    expectedZInPrevIm = cutOut3 + tempZShift;
    
    % Compute shift between two image stacks by cross-correlation
    crossCorr = [];
    for j = round(max(1,expectedZInPrevIm - 10)):round(min(prevImSize(3),expectedZInPrevIm + 10))
        corr2D = normxcorr2(overlapIm,prevIm(:,:,j));
        crossCorr = cat(3,crossCorr,corr2D);
    end
    % Find maximum peak of cross-correlation
    [~, idx] = max(crossCorr(:));
    [i1,i2,i3] = ind2sub(size(crossCorr),idx); % Subscripts from linear index
    % Convert peak position of cross-correlation to image shift values
    shift1 = i1 - size(overlapIm,1) + 1 - cutOutBegin1;
    shift2 = i2 - size(overlapIm,2) + 1 - cutOutBegin2;
    shift3 = i3 + max(1,expectedZInPrevIm - 10) - 1 - cutOut3;
    
    % Record accumulated shift values from first image stacks (apex end)
    shiftOfIms_pixel(i,:) = shiftOfIms_pixel(i-1,:) + [shift1 shift2 shift3];
    
    % Prepare for next iteration
    prevIm = currentIm;
    prevPoints = currentPoints;
    prevImSize = currentImSize;
end
shiftOfIms_pixel = round(shiftOfIms_pixel);

%% Merge overlapping regions among image stacks and adjust brightnesses

% Process image stacks one by one from the last stack (basal end)
load([SubFolderPath '\' TrimTif(imFileNames{end,1}) '.mat'],'-mat','originalIm');
processedIm = AdjustBrightness(originalIm, imBrightnessList(end)); % Brightness adjustment
save([SubFolderPath '\' TrimTif(imFileNames{end,1}) '.mat'],'-mat','processedIm','-append');

% Merge overlapping region between image stacks one by one. Image stack from
% apical side will be partially replaced with merged image. 
prevIm = processedIm;
prevImSize = size(prevIm);
for i = 2:size(imFileNames,1)
    j = size(imFileNames,1)-i+1;
    fprintf('Marging overlap between image stacks [%d and %d] ...\n',j,j+1)
    
    load([SubFolderPath '\' TrimTif(imFileNames{j,1}) '.mat'],'-mat','originalIm');
    currentIm = AdjustBrightness(originalIm, imBrightnessList(j));    
    currentImSize = size(currentIm);    
    shiftBtwn2 = shiftOfIms_pixel(j+1,:)-shiftOfIms_pixel(j,:); % Image shift
    
    % Cut out image of overlapped region from previous image
    prevBegin1 = max(1,1-shiftBtwn2(1,1));
    prevEnd1 = min(prevImSize(1),currentImSize(1)-shiftBtwn2(1,1));
    prevBegin2 = max(1,1-shiftBtwn2(1,2));
    prevEnd2 = min(prevImSize(2),currentImSize(2)-shiftBtwn2(1,2));
    prevBegin3 = max(1,1-shiftBtwn2(1,3));
    prevEnd3 = min(prevImSize(3),currentImSize(3)-shiftBtwn2(1,3));
    prevOverlapIm = prevIm(prevBegin1:prevEnd1,prevBegin2:prevEnd2,prevBegin3:prevEnd3);
    
    % Cut out image of overlapped region from current image
    currBegin1 = max(1,1+shiftBtwn2(1,1));
    currEnd1 = min(currentImSize(1),prevImSize(1)+shiftBtwn2(1,1));
    currBegin2 = max(1,1+shiftBtwn2(1,2));
    currEnd2 = min(currentImSize(2),prevImSize(2)+shiftBtwn2(1,2));
    currBegin3 = max(1,1+shiftBtwn2(1,3));
    currEnd3 = min(currentImSize(3),prevImSize(3)+shiftBtwn2(1,3));
    currentOverlapIm = currentIm(currBegin1:currEnd1,currBegin2:currEnd2,currBegin3:currEnd3);
    
    % Merge two 3D images
    % Compute plane with same distance from center coordinates of two image stacks
    overlapSize = size(currentOverlapIm);
    midpoint = overlapSize/2;
    % plane equation: ax + by + cz + d = 0
    a = shiftBtwn2(1,1); b = shiftBtwn2(1,2); c = shiftBtwn2(1,3);
    d = -a*midpoint(1)-b*midpoint(2)-c*midpoint(3);
    % Obtain all pixel positions of merging image
    tempIm = ones(overlapSize);
    xyzs = GetCoordOfPositivePixels(tempIm);
    % Compute distance between "xyz" and plane (ax + by + cz + d = 0);
    xyzDists = (a*xyzs(:,1)+b*(xyzs(:,2))+c*(xyzs(:,3))+d)/sqrt(a^2+b^2+c^2);
    xyzDists = (xyzDists - min(xyzDists(:))) / (max(xyzDists(:))-min(xyzDists(:))); % Normalize    
    % Compute contribution ratio for each image by logistic function with distance from the plane
    mixRatio = 1./(1+exp(-10*(xyzDists-0.5)));
    ratioMat = reshape(xyzDists,size(tempIm));
    margedIm = uint16(double(prevOverlapIm).*ratioMat + double(currentOverlapIm).*(1-ratioMat));
    
    % Replace overlap region in current image with marged image.
    processedIm = currentIm;
    processedIm(currBegin1:currEnd1,currBegin2:currEnd2,currBegin3:currEnd3) = margedIm;
    save([SubFolderPath '\' TrimTif(imFileNames{j,1}) '.mat'],'-mat','processedIm','-append');
    
    % Prepare for next iteration
    prevIm = processedIm;
    prevImSize = size(prevIm);
end

%% Combine curved lines along hair cells and vector list indicating inclinations
%at each point on the line.
%load lines and vectors of the first image stack
load([SubFolderPath '\' TrimTif(imFileNames{1,1}) '.mat'],'-mat','centerLinePerIm','inclinesPerIm');
centerLineOfCorti = centerLinePerIm;
inclineVectors = inclinesPerIm;

% Combine lines and vectors one by one
shiftOfIms = shiftOfIms_pixel*diag([PIXEL_WIDTH,PIXEL_WIDTH,Z_STEP]); % Convert to actual distance
for i = 2:size(imFileNames,1)
    load([SubFolderPath '\' TrimTif(imFileNames{i,1}) '.mat'],'-mat','centerLinePerIm','inclinesPerIm');
    % Compensate shift from the first image stack
    shiftedLine = centerLinePerIm + repmat(shiftOfIms(i,:),size(centerLinePerIm,1),1);
    [centerLineOfCorti, inclineVectors] ...
        = MergeLines(centerLineOfCorti, shiftedLine, inclineVectors, inclinesPerIm);
end

% Smooth line and vectors
inclineVectors = movmean(inclineVectors,25);
inclineVectors(:,3) = movmean(inclineVectors(:,3),20);

%% Create first linearized image
LINEAR_IM_WIDTH = 50; % Half size of first dimension in linearized image
LINEAR_IM_DEPTH = 25; % Half size of third dimension in linearized image
% Size of second dimension depends on cochlear sample

% Associate pixel positions in linearized image and Cartesian coordinates of original image stacks 
linearizedIm1 = ones(size(centerLineOfCorti,1),LINEAR_IM_WIDTH*2+1,LINEAR_IM_DEPTH*2+1);
xyzLinear = GetCoordOfPositivePixels(linearizedIm1); % Obtain pixel positions
xyzOrigin = ConvertLinearToOriginal(xyzLinear, centerLineOfCorti, inclineVectors ...
    , LINEAR_IM_WIDTH, LINEAR_IM_DEPTH);
xyzOrigin_pixel = xyzOrigin*diag([1/PIXEL_WIDTH,1/PIXEL_WIDTH,1/Z_STEP]); % Convert to pixel positions

% Convert pixel positions of linearized image to coordinates of original image stacks 
belongImList = zeros(size(xyzOrigin,1),size(imFileNames,1));
for i = 1:size(imFileNames,1)
    % Compensate shift from the first image stack
    xyzShifted = xyzOrigin_pixel - repmat(shiftOfIms_pixel(i,:),[size(xyzOrigin_pixel,1) 1]);
    isWithinIm = IsWithinIm(xyzShifted, imSizeList(i,:));
    belongImList(:,i) = (isWithinIm - sum(belongImList,2))>0; % Removing duplications
end

% Draw linearized image
linearizedIm1 = zeros(size(centerLineOfCorti,1),LINEAR_IM_WIDTH*2+1,LINEAR_IM_DEPTH*2+1);
for i = 1:size(imFileNames,1)
    fprintf('Making 1st linearized image [%d / %d] ...\n',i,numel(imFileNames))
    
    load([SubFolderPath '\' TrimTif(imFileNames{i,1}) '.mat'],'-mat','processedIm');
    % Obtain pixel values of linearized image from correspoinding image stacks
    linearizedIm1 = DrawLinearizedIm(linearizedIm1, processedIm, xyzLinear, xyzOrigin_pixel ...
        ,belongImList(:,i), imSizeList(i,:), shiftOfIms_pixel(i,:));
end
linearizedIm1 = flipud(permute(uint16(linearizedIm1),[2 1 3]));
ImWrite3D(linearizedIm1,[SubFolderPath '\linearizedIm1.tif']) % Export linearized image

%% Save data
if exist([SubFolderPath '\data.mat'],'file')
    save([SubFolderPath '\data.mat'],'linearizedIm1','counterClockWise','shiftOfIms_pixel' ...
        ,'PIXEL_WIDTH','Z_STEP','centerLineOfCorti','inclineVectors','imFileNames','imSizeList' ...
        ,'LINEAR_IM_WIDTH','LINEAR_IM_DEPTH','-append');
else
    save([SubFolderPath '\data.mat'],'linearizedIm1','counterClockWise','shiftOfIms_pixel' ...
        ,'PIXEL_WIDTH','Z_STEP','centerLineOfCorti','inclineVectors','imFileNames','imSizeList' ...
        ,'LINEAR_IM_WIDTH','LINEAR_IM_DEPTH')
end
disp('1st Linearization completed!')
