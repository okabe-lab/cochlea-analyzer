%*****************************************************************************
%This script is continuation of "main3.m". The script analyzes degree of 
%cell loss. The script get coordinates of inner and outer hair cells from 
%excel files "innerHairCells.xlsx" and "outerHairCells.xlsx", respectively. 
%These excel files are made by the preceding script "main3.m", and they can 
%be modified by manual inspection before running this script. This script uses
%a template data for 3D registration ("cochlearTemplate.mat"). 
%The result of anlysis will save in MAT-file "analyzeResults.mat". The script also
%makes "standardized_OuterHairCells.tif", a normalized image indicating locations 
%of outer hair cells.
%*****************************************************************************

if not(exist('mainPath','var'))
    mainPath = fileparts(mfilename('fullpath')); % Get fullpath of this script
end
cd(mainPath)

%% Information input
% Specify analysis range by x (horizontal) coordinates of linearized image
analysisRange = []; % for example [33 4860]
omitRange = []; % for example [3247 3330]

%% Get coordinates of inner hair cells
load([SubFolderPath '\data.mat'],'linearizedIm2','LINEAR_IM_WIDTH','LINEAR_IM_DEPTH','centerLineOfCorti2','inclineVectors2')
imSize = size(linearizedIm2);
tempxyz = xlsread(fullfile(SubFolderPath,'innerHairCells.xlsx'));
tempxyz = tempxyz(:,1:3);
innerHairCells = sortrows(SwitchColumn1_2(tempxyz),2);

%% Transform line along row of IHCs in linearized image into original image stack
innerStart = innerHairCells(1,:);
innerEnd = innerHairCells(end,:);
interval = 50;

startPointList = round(innerStart(1,2)+interval/2:interval:innerEnd(1,2))';
meanIHCs = zeros(numel(startPointList),3);
for i = 1:numel(startPointList)
    tempy = startPointList(i);
    f = (innerHairCells(:,2) > tempy-interval/2).*(innerHairCells(:,2) < tempy+interval/2);
    if sum(f)>0
        tempxyz = innerHairCells(f>0,:);
        meanIHCs(i,:) = mean(tempxyz);
    end
end
f = meanIHCs(:,1)==0;
meanIHCs(f>0,:) = [];
meanIHCs = [innerStart; meanIHCs; innerEnd];

% Compute spiral line along IHCs
meanIHCs2 = interp1(meanIHCs(:,2),meanIHCs,(innerStart(1,2):innerEnd(1,2))');
qline = [meanIHCs2(:,2) meanIHCs2(:,1) meanIHCs2(:,3)];
temp2 = qline - repmat([0,LINEAR_IM_WIDTH+1,LINEAR_IM_DEPTH+1],size(qline,1),1);
spiralLineAlongIHCs = centerLineOfCorti2(round(temp2(:,1)),:) - inclineVectors2(round(temp2(:,1)),:).*temp2(:,2);
spiralLineAlongIHCs(:,3) = spiralLineAlongIHCs(:,3) + temp2(:,3);

% Compute distance from apical end along spiral line for each IHC 
accumDistance = zeros(size(spiralLineAlongIHCs,1),1);
for j = 2:size(spiralLineAlongIHCs,1)
    accumDistance(j) = accumDistance(j-1) + norm(spiralLineAlongIHCs(j,:) - spiralLineAlongIHCs(j-1,:));
end
innerHairCellDistances = interp1(innerStart(1,2):innerEnd(1,2),accumDistance,innerHairCells(:,2),'linear','extrap');

%% Set cylindrical coordinate system along modiolus
interval = 50;
startPointList2 = (0:interval:accumDistance(end))';
pointsOnSpiral = zeros(size(startPointList2,1),3);
for i = 1:size(startPointList2,1)
    idx = knnsearch(accumDistance,startPointList2(i));
    pointsOnSpiral(i,:) = spiralLineAlongIHCs(idx,:);
end
% Compute initial value for modiolus estimation
center1 = ComputeCenterOfHelix([pointsOnSpiral startPointList2],0);
center2 = ComputeCenterOfHelix([pointsOnSpiral startPointList2],1);
centerDir = center1-center2;

vectorA = cross(centerDir/norm(centerDir),[0,0,1]);
deg = acos(dot(centerDir/norm(centerDir),[0,0,1]));
rotatedSpiral1 = RodriguesRotation(pointsOnSpiral(:,1:3),vectorA,deg);

% Estimation of modioulus
f = @(x) ComputeRssOfCochlearHelixFit(rotatedSpiral1,x(1),x(2));
rotateAngles = fminsearch(f,[0,0]); % Search best angles for fitting

% Transform IHCs' coordinates to new coordinate system 
temp = RodriguesRotation(rotatedSpiral1,[1,0,0],rotateAngles(1));
rotatedSpiral2 = RodriguesRotation(temp,[0,1,0],rotateAngles(2));
[~,cent1,cent2] =ComputeRssOfCochlearHelixFit(rotatedSpiral1,rotateAngles(1),rotateAngles(2));

temp2 = round(SwitchColumn1_2(innerHairCells)) - repmat([0,LINEAR_IM_WIDTH+1,LINEAR_IM_DEPTH+1],size(innerHairCells,1),1);
qline = centerLineOfCorti2(temp2(:,1),:) - inclineVectors2(temp2(:,1),:).*temp2(:,2);
qline(:,3) = qline(:,3) + temp2(:,3);
zList2 = RodriguesRotation(qline,vectorA,deg);
temp2 = RodriguesRotation(zList2,[1,0,0],rotateAngles(1));
temp3 = RodriguesRotation(temp2,[0,1,0],rotateAngles(2));
temp4 = [1 0 0 -cent1; 0 1 0 -cent2; 0 0 1 -rotatedSpiral2(1,3); 0 0 0 1]*[temp3 ones(size(temp3,1),1)]';
originalIHCs1 = temp4(1:3,:)';
if originalIHCs1(end,3)<0 % Adjust Z direction
    originalIHCs1(:,3)=-originalIHCs1(:,3);
    originalIHCs1(:,1)=-originalIHCs1(:,1);
end

% Set azimuth = 0 at 1st IHC and calculate azimuth for each IHC
azimuthIHC1 = zeros(size(originalIHCs1,1),1);
radialDistIHC1 = zeros(size(originalIHCs1,1),1);
for i = 1:size(originalIHCs1,1)
    temp = originalIHCs1(i,1:2);
    radialDistIHC1(i,1) = norm(temp);    
    temp = temp/norm(temp);
    if i == 1
        v0 = temp;
        v00 = temp;
        azimuthIHC1(i,1) = 0;
        judge = cross([v00 0],[temp 0]);
        continue
    end
    
    temp2 = acos(dot(v0,temp));    
    v0 = temp;
    azimuthIHC1(i,1) = azimuthIHC1(i-1,1) + temp2;
end
originalIHCs2 = [radialDistIHC1.*cos(azimuthIHC1) radialDistIHC1.*sin(azimuthIHC1) originalIHCs1(:,3)];

%% 3D Registration with template data
load('cochlearTemplate.mat','z_temp','s_temp')
azimuthList0 = (0:0.01:pi*4)';
midAzimuth = max(s_temp)/2;

temp = interp1(azimuthIHC1,originalIHCs2(:,3),azimuthList0);
temp(isnan(temp),:) = [];
zList = temp;

templateZList = z_temp(100:min(numel(zList),numel(z_temp))-50,:);
r = normxcorr2(templateZList,zList); % Calculate lag of ƒÓ(z)
[~,I] = max(r);
lag = I - numel(templateZList) - 99;
azimuthIHC2 = azimuthIHC1 - lag*0.01;

% Obtain intersection points
zList2 = NaN(numel(s_temp),1);
st1 = max(1,lag);
ed1 = min(numel(s_temp)+lag-1,numel(zList));
st2 = max(1,-lag);
zList2(st2:st2+(ed1-st1)) = zList(st1:ed1);

% Adjust z coordinates
temp = z_temp-zList2;
f = @(x) nansum((temp - x).^2);
zlag = fminsearch(f,0);
adjustedZ = originalIHCs2(:,3) + zlag;
midOfSpiral = interp1(azimuthIHC2,innerHairCellDistances,midAzimuth);

%% Get coordinates of outer hair cells
xyz = xlsread([SubFolderPath '\outerHairCells.xlsx']);
outerHairCells = sortrows(SwitchColumn1_2(xyz),2);
outerHairCellDistances = interp1(innerStart(1,2):innerEnd(1,2),accumDistance,outerHairCells(:,2),'linear','extrap');

%% Analyze cell loss
% Set x range for analysis
if ~isempty(omitRange)
    f = zeros(size(outerHairCellDistances));
    for k = 1:size(omitRange,1)
        com = interp1(innerHairCells(:,2),innerHairCellDistances,omitRange(k,:));
        f((outerHairCellDistances>com(1)).*(outerHairCellDistances<=com(2))>0)=1;
    end
    outerHairCellDistances(f>0) = [];
end

if isempty(analysisRange)
    startPoint = round(min(innerHairCells(:,2)));
    endPoint = round(max(outerHairCells(:,2)));
else
    startPoint = min(analysisRange);
    endPoint = max(analysisRange);
end

% Cell loss detection
[standardizedIm, simplane2, spacePixelList, spaceCenterList, S] ...
    = DetectCellLoss(outerHairCells, startPoint, endPoint, imSize, omitRange);
emptylist = cat(1, S.Area);
lossClusterSizes = round(emptylist/25);

%% Compute adjusted distances from basal end
distanceAdjustment = 3200;
distanceIHCs = -(innerHairCellDistances - midOfSpiral)+distanceAdjustment;
distanceOHCs = -(outerHairCellDistances - midOfSpiral)+distanceAdjustment;
distanceLostPixels = interp1(innerHairCells(:,2),distanceIHCs,spacePixelList(:,2));
endOfRange = interp1(innerHairCells(:,2),distanceIHCs,endPoint,'linear','extrap');
distanceLostCenters = interp1(innerHairCells(:,2),distanceIHCs,spaceCenterList(:,2));

if ~isempty(omitRange)
    com = omitRange;
    fi = zeros(size(distanceIHCs));
    fe = zeros(size(distanceLostPixels));
    fo = zeros(size(distanceOHCs));
    for k = 1:size(omitRange,1)
        com(k,:) = interp1(innerHairCells(:,2),distanceIHCs,omitRange(k,[2 1]));
        fi((distanceIHCs>com(k,1)).*(distanceIHCs<=com(k,2))>0)=1;
        fe((distanceLostPixels>com(k,1)).*(distanceLostPixels<=com(k,2))>0)=1;
        fo((distanceOHCs>com(k,1)).*(distanceOHCs<=com(k,2))>0)=1;
    end
    distanceIHCs(fi>0) = [];
    distanceLostPixels(fe>0) = [];
    distanceOHCs(fo>0) = [];
else
    com = [];
end

distanceIHCs(distanceIHCs<endOfRange)=[];
distanceOHCs(distanceOHCs<endOfRange)=[];
distanceLostPixels(distanceLostPixels<endOfRange)=[];
lossClusterSizes = lossClusterSizes(distanceLostCenters>=endOfRange);

%% Count number of cell loss in segments 
interval = 50;
startPointList3 = 0:interval:6000;

lostCellNo = NaN(numel(startPointList3)-1,1);
for i = 1:numel(lostCellNo)
    tempst = startPointList3(i);
    temped = startPointList3(i+1);
    % Check number of cells in segment
    f1 = distanceIHCs >= tempst;
    f2 = distanceIHCs < temped;
    tin = distanceIHCs(f1.*f2>0);
    f3 = distanceOHCs >= tempst;
    f4 = distanceOHCs < temped;
    tout = distanceOHCs(f3.*f4>0);
    xCoords = [tin; tout];
    if isempty(xCoords)
        continue
    elseif (max(xCoords)-min(xCoords)) < interval*0.5
        continue
    end
    
    % Count cell loss
    f5 = distanceLostPixels >= tempst;
    f6 = distanceLostPixels < temped;
    emptyPixelNo = distanceLostPixels(f5.*f6>0);
    lostCellNo(i) = round(numel(emptyPixelNo)/25);
end
lostCellNoList50 = lostCellNo;

%% Calculate total loss ratio 
conum = numel(distanceOHCs);
totalLossRatio = sum(lossClusterSizes)/(conum+sum(lossClusterSizes));

%% Save and export results
tempim = cat(3,standardizedIm,simplane2,zeros(size(standardizedIm)));
imwrite(tempim(2:16,:,:),[SubFolderPath '\standardized_OuterHairCells.tif'])
save([SubFolderPath '\analyzeResults.mat'],'lostCellNoList50','totalLossRatio','lossClusterSizes');

disp('main4.m completed!')