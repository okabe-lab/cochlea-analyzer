%*****************************************************************************
%This script is continuation of "main1.m". This script perform detection of
%inner hair cells on 1st linearized image, and create 2nd linearized image 
%based on the detected inner hair cell row ("linearizedIm2.tif"). It also 
%create an image with estimated inner hair cell locations ("innerHairCells.tif"), 
%and excel file with their coordinates("innerHairCells.xlsx"). This script uses 
%machine learning models ("machineLearningModels.mat").  
%*****************************************************************************

clearvars -except SubFolderPath
if not(exist('mainPath','var'))
    mainPath = fileparts(mfilename('fullpath')); % Get fullpath of this script
end
cd(mainPath)
load([SubFolderPath '\data.mat'],'linearizedIm1','shiftOfIms_pixel','PIXEL_WIDTH','Z_STEP' ...
    ,'centerLineOfCorti','inclineVectors','imFileNames','imSizeList','LINEAR_IM_WIDTH','LINEAR_IM_DEPTH')
load('machineLearningModels.mat','imdl1s','imdl2');
imSize = size(linearizedIm1);

%% Detection of inner hair cells from 1st image
coordIHCs1 = DetectInnerHairCells(linearizedIm1,imdl1s{1,1},imdl2);
coordIHCs1(:,3) = movmean(coordIHCs1(:,3),4);
if coordIHCs1(end,2) < imSize(2)-45
    coordIHCs1 = [coordIHCs1; coordIHCs1(end,1) imSize(2) coordIHCs1(end,3)];
end
disp('Inner hair cells in 1st image detected...')

%% Create 2nd linearized image
% Transform line along IHCs of linearized image to spiral in original image
tempLine = interp1(coordIHCs1(:,2),[coordIHCs1(:,1) coordIHCs1(:,3)],round(min(coordIHCs1(:,2)):max(coordIHCs1(:,2))) ...
    ,'linear','extrap');
tempLine = [round(min(coordIHCs1(:,2)):max(coordIHCs1(:,2)))' tempLine];
tempLine2 = tempLine - repmat([0,LINEAR_IM_WIDTH+1,LINEAR_IM_DEPTH+1],size(tempLine,1),1);
spiralAlongIHCs = centerLineOfCorti(tempLine2(:,1),:) - inclineVectors(tempLine2(:,1),:).*tempLine2(:,2);
spiralAlongIHCs(:,3) = spiralAlongIHCs(:,3) + tempLine2(:,3);

% Measure distance along spiral
accumDistance = zeros(size(spiralAlongIHCs,1),1);
for j = 2:size(spiralAlongIHCs,1)
    accumDistance(j) = accumDistance(j-1) + norm(spiralAlongIHCs(j,:) - spiralAlongIHCs(j-1,:));
end

% Obtain points along IHCs' line/spiral with 1 micrometer intervals
interval = 40;
startPointList = (0:interval:accumDistance(end)-interval)';
pointsAlongIHCLine = zeros(size(startPointList,1),3);
for j = 1:size(startPointList,1)
    idx = knnsearch(accumDistance,startPointList(j));
    pointsAlongIHCLine(j,:) = tempLine(idx,:); 
end
pointsAlongIHCLine = SwitchColumn1_2(pointsAlongIHCLine);
pointsAlongIHCLine2 = pointsAlongIHCLine;
pointsAlongIHCLine2(:,1) = pointsAlongIHCLine(:,1)+30; %Shift toward center of linearized image

% Adjust Z coordinates of points
width = 40;
interval = 6;
for j = 1:size(pointsAlongIHCLine,1)
    temp = pointsAlongIHCLine(j,:);
    pointsAlongIHCLine2(j,3) = AdjustZCoord(linearizedIm1,temp,width,interval);       
end
pointsAlongIHCLine2 = (pointsAlongIHCLine+pointsAlongIHCLine2)/2;

centerLine1 = interp1(pointsAlongIHCLine2(:,2),[pointsAlongIHCLine2(:,1) pointsAlongIHCLine2(:,3)],1:size(linearizedIm1,2) ...
    ,'linear','extrap');
centerLine1 = [(1:size(linearizedIm1,2))' centerLine1];
tempLine2 = centerLine1 - repmat([0,LINEAR_IM_WIDTH+1,LINEAR_IM_DEPTH+1],size(centerLine1,1),1);
centerSpiral1 = centerLineOfCorti(tempLine2(:,1),:) - inclineVectors(tempLine2(:,1),:).*tempLine2(:,2);
centerSpiral1(:,3) = centerSpiral1(:,3) + tempLine2(:,3);

% Measure distance along center spiral
accumDistance = zeros(size(centerSpiral1,1),1);
for j = 2:size(centerSpiral1,1)
    accumDistance(j) = accumDistance(j-1) + norm(centerSpiral1(j,:) - centerSpiral1(j-1,:));
end

% Obtain points along center line/spiral of Corti in linearized image with 1 micrometer intervals
interval = 1;
startPointList = (0:interval:accumDistance(end))';
centerLineOfCorti2 = zeros(size(startPointList,1),3);
centerLine2 = zeros(size(startPointList,1),3);
for j = 1:size(startPointList,1)
    idx = find((accumDistance-startPointList(j))>0,1);
    tempdif = accumDistance(idx)-startPointList(j);
    if tempdif > 0
        centerLineOfCorti2(j,:) = centerSpiral1(idx,:) + (centerSpiral1(idx-1,:)-centerSpiral1(idx,:))*tempdif;
    else 
        centerLineOfCorti2(j,:) = centerSpiral1(end,:);
    end
    if tempdif > 0
        tempvect = (centerLine1(idx-1,:)-centerLine1(idx,:));
        tempvect = (tempvect/norm(tempvect))*tempdif;
        centerLine2(j,:) = centerLine1(idx,:) + tempvect;
    else 
        centerLine2(j,:) = centerLine1(end,:);
    end
end

% Re-obtain points along IHC line/spiral of Corti in linearized image with 1 micrometer intervals
pointsAlongIHCLine3 = interp1(pointsAlongIHCLine(:,2),[pointsAlongIHCLine(:,1) pointsAlongIHCLine(:,3)],centerLine2(:,1) ...
    ,'linear','extrap');
pointsAlongIHCLine3 = [centerLine2(:,1) pointsAlongIHCLine3];
temp = centerLine2(1,:);
tempLine2 = pointsAlongIHCLine3(1:50,:)-temp;
temp3 = centerLine2(1,:)-centerLine2(50,:);
if temp3(2)<0 % Avoid erroneous detection around apical end
    [~,idx] = min(abs(dot(repmat(temp3,size(tempLine2,1),1),tempLine2,2)));
    tinline_3 = pointsAlongIHCLine3(idx:50,:);
    tinline_3 = interp1([1;50],[tinline_3(1,:); tinline_3(end,:)],1:50);
    pointsAlongIHCLine3(1:50,:) = tinline_3;
end
tempLine2 = pointsAlongIHCLine3 - repmat([0,LINEAR_IM_WIDTH+1,LINEAR_IM_DEPTH+1],size(pointsAlongIHCLine3,1),1);
temp3 = interp1(1:size(centerLineOfCorti,1),centerLineOfCorti,pointsAlongIHCLine3(:,1),'linear','extrap');
temp4 = interp1(1:size(inclineVectors,1),inclineVectors,pointsAlongIHCLine3(:,1),'linear','extrap');
spiralAlongIHCs = temp3 - temp4.*tempLine2(:,2);
spiralAlongIHCs(:,3) = spiralAlongIHCs(:,3) + tempLine2(:,3);

% Adjust distances among points on spiral around apical ends
num = 300;
accumDistance = zeros(num,1);
for j = 2:num
    accumDistance(j) = accumDistance(j-1) + norm(spiralAlongIHCs(j,:) - spiralAlongIHCs(j-1,:));
end
interval = accumDistance(end)/num;
startPointList = (0:interval:accumDistance(end))';
tempSpiral = zeros(num,3);
for j = 1:num
    idx = find((accumDistance-startPointList(j))>0,1);
    tempdif = accumDistance(idx)-startPointList(j);
    if tempdif > 0
        tempSpiral(j,:) = spiralAlongIHCs(idx,:) + (spiralAlongIHCs(idx-1,:)-spiralAlongIHCs(idx,:))*tempdif;
    else 
        tempSpiral(j,:) = spiralAlongIHCs(num,:);
    end
end
spiralAlongIHCs(1:num,:) = tempSpiral;

num = 300;
accumDistance = zeros(num,1);
for j = 2:num
    accumDistance(j) = accumDistance(j-1) + norm(spiralAlongIHCs(num+j,:) - spiralAlongIHCs(num+j-1,:));
end
interval = accumDistance(end)/num;
startPointList = (0:interval:accumDistance(end))';
tempSpiral = zeros(num,3);
for j = 1:num
    idx = find((accumDistance-startPointList(j))>0,1);
    tempdif = accumDistance(idx)-startPointList(j);
    if tempdif > 0
        tempSpiral(j,:) = spiralAlongIHCs(num+idx,:) + (spiralAlongIHCs(num+idx-1,:)-spiralAlongIHCs(num+idx,:)) ...
            *tempdif;
    else 
        tempSpiral(j,:) = spiralAlongIHCs(num,:);
    end
end
spiralAlongIHCs(num+1:num*2,:) = tempSpiral;

% Compute inclinations from two spirals
vector2 = spiralAlongIHCs-centerLineOfCorti2;
vector2 = vector2./sqrt(sum(vector2.^2,2));
inclineVectors2 = movmean(vector2,7);

%% Create linearized image based on center spiral line along organ of Corti and inclinations
Z_STEP2 = 1;
linearizedIm2 = ones(size(centerLineOfCorti2,1),LINEAR_IM_WIDTH*2+1,LINEAR_IM_DEPTH*2+1);
xyzLinear = GetCoordOfPositivePixels(linearizedIm2);
xyzOrigin = ConvertLinearToOriginal(xyzLinear, centerLineOfCorti2, inclineVectors2 ...
    , LINEAR_IM_WIDTH, LINEAR_IM_DEPTH);
xyzOrigin_pixel = xyzOrigin*diag([1/PIXEL_WIDTH,1/PIXEL_WIDTH,1/Z_STEP]); % Convert to pixel positions

belongImList = zeros(size(xyzOrigin,1),size(imFileNames,1));
for i = 1:size(imFileNames,1)
    % Compensate shift from the first image stack
    xyzShifted = xyzOrigin_pixel - repmat(shiftOfIms_pixel(i,:),[size(xyzOrigin_pixel,1) 1]);
    isWithinIm = IsWithinIm(xyzShifted, imSizeList(i,:));
    belongImList(:,i) = (isWithinIm - sum(belongImList,2))>0; % Removing duplications
end

linearizedIm2 = zeros(size(centerLineOfCorti2,1),LINEAR_IM_WIDTH*2+1,LINEAR_IM_DEPTH*2+1);
for j = 1:size(imFileNames,1)
    fprintf('Making 2nd linearized image [%d / %d] ...\n',j,numel(imFileNames))
    load([SubFolderPath '\' TrimTif(imFileNames{j,1}) '.mat'],'-mat','processedIm');
    linearizedIm2 = DrawLinearizedIm(linearizedIm2,processedIm,xyzLinear,xyzOrigin_pixel,belongImList(:,j),imSizeList(j,:) ...
        ,shiftOfIms_pixel(j,:));
end
linearizedIm2 = flipud(permute(uint16(linearizedIm2),[2 1 3]));

ImWrite3D(linearizedIm2,[SubFolderPath '\linearizedIm2.tif']);
disp('Re_linearization completed...')

%% Detection of inner hair cells in 2nd image
innerCells2 = DetectInnerHairCells(linearizedIm2,imdl1s{1,1},imdl2);
innerCells2(:,3) = movmean(innerCells2(:,3),4);

markedImOfIHC = DrawMarkedIm(linearizedIm2,round(innerCells2),2,1);
markedImOfIHC = uint16(markedImOfIHC*1000);
disp('Inner hair cells in 2nd image detected!')

%% Save and export data
save([SubFolderPath '\data.mat'],'linearizedIm2','centerLineOfCorti2','inclineVectors2','innerCells2','-append')
xlswrite([SubFolderPath '\innerHairCells.xlsx'],sortrows(SwitchColumn1_2(round(innerCells2))),1)
ImWrite3D(markedImOfIHC,[SubFolderPath '\innerHairCells.tif']);
