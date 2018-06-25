%*****************************************************************************
%This script is continuation of "main2.m". This script perform detection of
%outer hair cells on 2nd linearized image. It creat an image with estimated
%outer hair cell locations ("outerHairCells.tif"), excel files with estimated
%locations of inner and outer hair cells ("outerHairCells.xlsx"). This script uses 
%machine learning models ("machineLearningModels.mat").
%*****************************************************************************

clearvars -except SubFolderPath
if not(exist('mainPath','var'))
    mainPath = fileparts(mfilename('fullpath')); % Get fullpath of this script
end
cd(mainPath)
load([SubFolderPath '\data.mat'],'linearizedIm2','innerCells2')
load('machineLearningModels.mat','omdl1s','omdl2','net1','net2');

%% First selection of candidate locations of outer hair cells
Mdl = omdl1s{1,1}; % Remove noise
[predictor1,~,predictionScores1,peakGroupStats, labeledIm, correlationMatrix] = PredictOuterHairCells(linearizedIm2,innerCells2,Mdl);
disp('Cell candidates obtained...')

imSize = size(linearizedIm2);
width = 200;
num1 = round(imSize(2)/width);
num2 = round(width/2);
startPointList = [0:floor(imSize(2)/num1):floor(imSize(2)/num1)*(num1-1) imSize(2)];
idxCell = cell(num1,1);
for j = 1:num1
    f = (predictor1(:,6) >= startPointList(j)).*(predictor1(:,6) < startPointList(j+1));
    tempScore = predictionScores1(:,2).*f;
    [~,tempIdx] = sort(tempScore,'descend');
    idxCell{j,1} = tempIdx(1:num2,:);
end
selectedIndexList = cat(1,idxCell{:});

%% Additional feature quantity extraction
subPredict1 = predictor1(selectedIndexList,:);
selectedPredScore1 = predictionScores1(selectedIndexList,2);

% Search neighbor "cell candidates"
xyzList1 = subPredict1(:,5:7);
selectedPixelList = subPredict1(:,8:10);
xyzList3 = subPredict1(:,11:13);
[distanceList,xDistList,yDistList,zDistList] = ComputeMinDistBtwnClusters(xyzList1,selectedPixelList,xyzList3);

% Add neighbors data
subPredict2 = zeros(size(subPredict1,1),24);
neighborIndex = zeros(size(subPredict1,1),7);
isNeighbor = (distanceList > 4).*(distanceList <16);
for j = 1:size(subPredict1,1)
    tempNeighbor = find(isNeighbor(:,j));
    tempVectorList = [xDistList(tempNeighbor,j), yDistList(tempNeighbor,j), zDistList(tempNeighbor,j)];
    normal = sqrt(sum(tempVectorList.^2,2));
    normalizedVectList = tempVectorList./normal;
    angleList = acos(normalizedVectList(:,1));
    f = normalizedVectList(:,2)<0;
    angleList(f>0) = 2*pi - angleList(f>0);
    tempScore = selectedPredScore1(tempNeighbor);
    
    for k = 1:6
        f1 = angleList >= (k-1)*(pi/3);
        f2 = angleList < k*(pi/3);
        [maxValue,idx] = max(tempScore.*f1.*f2);
        if maxValue >0
            subPredict2(j,(k-1)*4+1:k*4)=[tempVectorList(idx,:),tempScore(idx)];
            neighborIndex(j,k) = tempNeighbor(idx);
        end
    end
end
[sortedDistance,I] = sort(distanceList,2);
neighborIndex2 = I(:,2);
distanceList2 = sortedDistance(:,2);
num = size(distanceList,1);
tempIdx = sub2ind([num num],(1:num)',neighborIndex2);
vectorList = [xDistList(tempIdx) yDistList(tempIdx), zDistList(tempIdx)];
subPredict3 = [vectorList selectedPredScore1(neighborIndex2) distanceList2];
neighborIndex(:,7) = neighborIndex2;

% Add image around each cell candidate
expandedIm = padarray(linearizedIm2,[100 100 5]);
subPredict4 = zeros(size(subPredict1,1),21*7);
for j = 1:size(subPredict1,1)
    temp = round(subPredict1(j,5:7));
    stx = temp(1)-28+100; edx = stx+63-1;
    sty = temp(2)-10+100; edy = sty+21-1;
    stz = temp(3)-1+5; edz = stz+3-1;
    tempIm = histeq(max(expandedIm(stx:edx,sty:edy,stz:edz),[],3));
    tempIm = imresize(tempIm,[21,7]);
    tempIm = double(tempIm);
    subPredict4(j,:) = tempIm(:)/prctile(tempIm(:),90);
end

predictor2 = [subPredict1 subPredict2 subPredict3 subPredict4];

% Obtain coordinates of cell candidates
cellCandidCoord = predictor2(:,5:7);
for j = 1:size(selectedIndexList)
    temp = SwitchColumn1_2(peakGroupStats(selectedIndexList(j)).PixelList);
    [idx,~] = knnsearch(temp,cellCandidCoord(j,:)); % Search nearest pixel
    cellCandidCoord(j,:) = temp(idx,:);
end

% Image data for "net1"
wx = 69; wy = 39; wz = 3;
predictorImList = zeros(69,39,1,size(cellCandidCoord,1));
for j = 1:size(cellCandidCoord,1)
    temp = round(cellCandidCoord(j,:));
    stx = temp(1)-floor(wx/2)+100; edx = stx+wx-1;
    sty = temp(2)-floor(wy/2)+100; edy = sty+wy-1;
    stz = temp(3)-floor(wz/2)+5; edz = stz+wz-1;
    tempIm = double(max(expandedIm(stx:edx,sty:edy,stz:edz),[],3));
    tempIm = histeq(tempIm/max(tempIm(:)));
    predictorImList(:,:,:,j) = tempIm;
end

%% Second selection of candidate locations of inner hair cells
Mdl = omdl2;
[binaryPredictions1,predictionScores2] = predict(Mdl,predictor2);
rowPredictions = double(classify(net1,predictorImList));     
selectedPixelList = SwitchColumn1_2(cat(1,peakGroupStats(selectedIndexList).PixelList));
imSize = size(linearizedIm2);
rowPredictions(binaryPredictions1==-1)=-1;
cellCandidData = [selectedIndexList rowPredictions cellCandidCoord predictionScores2(:,2)];

leftEnd = round(min(innerCells2(:,2)));
f1 = cellCandidData(:,4) < (leftEnd+200);
f2 = binaryPredictions1 == 1;
apicalCandidates = cellCandidData(f1.*f2>0,:);

%% Detect rows of outer hair cells
disp('Estimating outer hair cells...')
% Remove candidates near boundaries
for j = 1:size(cellCandidData,1)
    tempData = cellCandidData(j,:);
    if tempData(2) ~=-1
        if tempData(4) > imSize(2)-10 || tempData(5) < 6 || tempData(5) > imSize(1)-6
            cellCandidData(j,2) = -1;
        else
            nearestPeakPoint = round(tempData(3:5));
            if nearestPeakPoint(1)+5>imSize(1)
                cellCandidData(j,2) = -1;
            else              
                pixelValue1 = linearizedIm2(nearestPeakPoint(1),nearestPeakPoint(2)+10,nearestPeakPoint(3));
                pixelValue2 = linearizedIm2(nearestPeakPoint(1)-5,nearestPeakPoint(2),nearestPeakPoint(3));
                pixelValue3 = linearizedIm2(nearestPeakPoint(1)+5,nearestPeakPoint(2),nearestPeakPoint(3));
            if pixelValue1 == 0 || pixelValue2 == 0 || pixelValue3 == 0
                cellCandidData(j,2) = -1;
            end
            end
        end
    end
end

% Remove too close candidate pairs
for rowNo = 1:3
    row1List = sortrows(cellCandidData(cellCandidData(:,2)==1,3:5),2);
    row2List = sortrows(cellCandidData(cellCandidData(:,2)==2,3:5),2);
    row3List = sortrows(cellCandidData(cellCandidData(:,2)==3,3:5),2);
    allRowList = {row1List;row2List;row3List};
    row1Intervals = diff(row1List);
    row2Intervals = diff(row2List);
    row3Intervals = diff(row3List);
    allRowIntervals = {row1Intervals;row2Intervals;row3Intervals};
    
    for j = 1:size(allRowIntervals{rowNo,1})
        if allRowIntervals{rowNo,1}(j,2) < 5
            xyz1 = allRowList{rowNo,1}(j,:);
            xyz2 = allRowList{rowNo,1}(j+1,:);
            f1 = (cellCandidCoord(:,1) == xyz1(1)) .* (cellCandidCoord(:,2) == xyz1(2)) .* (cellCandidCoord(:,3) == xyz1(3));
            f2 = (cellCandidCoord(:,1) == xyz2(1)) .* (cellCandidCoord(:,2) == xyz2(2)) .* (cellCandidCoord(:,3) == xyz2(3));
            if abs(xyz1(1)-xyz2(1))<4
                cellCandidData(f1>0,2) = -1;
                cellCandidData(f2>0,2) = -1;
            else
                cellCandidData(f1>0,2) = 4;
                cellCandidData(f2>0,2) = 4;
            end
        end
    end
end

% Recover removed candidates with high prediction scores
f1 = cellCandidData(:,2)==-1;
f2 = cellCandidData(:,6)>0.8;
cellCandidData(f1.*f2>0,2) = 4;

% Recover removed candidates in gap space
for rowNo = 1:3
    row1List = sortrows(cellCandidData(cellCandidData(:,2)==1,3:5),2);
    row2List = sortrows(cellCandidData(cellCandidData(:,2)==2,3:5),2);
    row3List = sortrows(cellCandidData(cellCandidData(:,2)==3,3:5),2);
    allRowList = {row1List;row2List;row3List};
    row1Intervals = diff(row1List);
    row2Intervals = diff(row2List);
    row3Intervals = diff(row3List);
    allRowIntervals = {row1Intervals;row2Intervals;row3Intervals};
    cellCandidCoord = cellCandidData(:,3:5);
    isRemovedList = cellCandidData(:,2) ~= -1;
    
    for j = 1:size(allRowIntervals{rowNo,1})
        if allRowIntervals{rowNo,1}(j,2) > 12 && allRowIntervals{rowNo,1}(j,2) < 50
            xyz1 = allRowList{rowNo,1}(j,:);
            xyz2 = allRowList{rowNo,1}(j+1,:);
            gapCoords = interp1([xyz1(2);xyz2(2)],[xyz1;xyz2],ceil(xyz1(2)+4):floor(xyz2(2)-4));
            
            f1 = (cellCandidCoord(:,2) > allRowList{rowNo,1}(j,2)).*(cellCandidCoord(:,2) < allRowList{rowNo,1}(j+1,2));
            f2 = (cellCandidCoord(:,1) > (min(allRowList{rowNo,1}(j,1),allRowList{rowNo,1}(j+1,1))-8)) ...
                .* (cellCandidCoord(:,1) < (max(allRowList{rowNo,1}(j,1),allRowList{rowNo,1}(j+1,1))+8));
            recoverCandidIndList = find(isRemovedList.*f1.*f2);
            
            if ~isempty(recoverCandidIndList)
                idxList = cellCandidData(recoverCandidIndList,1);
                candidPixelList = SwitchColumn1_2(cat(1,peakGroupStats(idxList).PixelList));
                [idxs3,nearestDist] = knnsearch(candidPixelList,gapCoords);
                candidPixelList2 = candidPixelList(idxs3(nearestDist<3,:),:);
                if ~isempty(candidPixelList2)
                    linearIndice = sub2ind(imSize,candidPixelList2(:,1),candidPixelList2(:,2),candidPixelList2(:,3));
                    groupNo1 = labeledIm(linearIndice);
                    groupNo2 = unique(groupNo1);
                    for kk = size(groupNo2,1)
                        cellCandidData(cellCandidData(:,1) == groupNo2(kk),2) = rowNo;
                        isGroupList = groupNo1 == groupNo2(kk);
                        tempCoords = candidPixelList2(isGroupList>0,:);                        
                        cellCandidData(cellCandidData(:,1) == groupNo2(kk),3:5) = tempCoords(1,:);
                    end
                end
            end
        end
    end
end

rightEnd = round(max(innerCells2(:,2)));
cellCandidCoord = cellCandidData(:,3:5);

% Compute average z coordinate in segments
interval = 60;
width = 60;
startPointList = [(leftEnd:interval:rightEnd-width), rightEnd-width+1]';
zList = zeros(size(startPointList,1),4);
for j = 1:size(startPointList,1)
    f1 = (cellCandidData(:,2)==1).*(cellCandidCoord(:,2) >= startPointList(j)).*(cellCandidCoord(:,2) < startPointList(j)+width);
    f2 = (cellCandidData(:,2)==2).*(cellCandidCoord(:,2) >= startPointList(j)).*(cellCandidCoord(:,2) < startPointList(j)+width);
    f3 = (cellCandidData(:,2)==3).*(cellCandidCoord(:,2) >= startPointList(j)).*(cellCandidCoord(:,2) < startPointList(j)+width);
    row1Zs = cellCandidCoord(f1>0,3);
    row2Zs = cellCandidCoord(f2>0,3);
    row3Zs = cellCandidCoord(f3>0,3);
    zList(j,:) = [startPointList(j)+width/2 mean(row1Zs) mean(row2Zs) mean(row3Zs)];
end
f = sum(isnan(zList),2)>0;
zList(f>0,:) = [];
interpZList = interp1(zList(:,1),zList(:,2:4),(1:imSize(2))','linear','extrap');

% Remove candidates far from average value of z coordinate
for j = 1:size(cellCandidData,1)
    if cellCandidData(j,2)==1 || cellCandidData(j,2)==2 || cellCandidData(j,2)==3
        if abs(cellCandidData(j,5)-interpZList(round(cellCandidData(j,4)),cellCandidData(j,2)))>9
            cellCandidData(j,2) = -1;
        end
    end
end

% Compute average inner-row intervals 
interval = 250;
width = 250;
startPointList = [(leftEnd:interval:rightEnd-width), rightEnd-width+1]';
innerRowIntervals = zeros(size(startPointList,1),2);
for j = 1:size(startPointList,1)
    f1 = (cellCandidData(:,2)==1).*(cellCandidCoord(:,2) >= startPointList(j)).*(cellCandidCoord(:,2) < startPointList(j)+width);
    f2 = (cellCandidData(:,2)==2).*(cellCandidCoord(:,2) >= startPointList(j)).*(cellCandidCoord(:,2) < startPointList(j)+width);
    f3 = (cellCandidData(:,2)==3).*(cellCandidCoord(:,2) >= startPointList(j)).*(cellCandidCoord(:,2) < startPointList(j)+width);
    temp1 = diff(sort(cellCandidCoord(f1>0,2)));
    temp2 = diff(sort(cellCandidCoord(f2>0,2)));
    temp3 = diff(sort(cellCandidCoord(f3>0,2)));
    innerRowIntervals(j,:) = [startPointList(j)+width/2 median([temp1; temp2; temp3])];
end
f = sum(isnan(innerRowIntervals),2)>0;
innerRowIntervals(f>0,:) = [];
innerRowIntervals(:,2) = medfilt1(innerRowIntervals(:,2));
interpIntervalList = interp1(innerRowIntervals(:,1),innerRowIntervals(:,2),(1:imSize(2))','linear','extrap');

% Adjust z coordinates of cell candidates
for j = 1:size(cellCandidData,1)
    if cellCandidData(j,2)~=-1
        tempRowLabel = cellCandidData(j,2);
        if tempRowLabel == 4
            tempRowLabel = 3;
        end
        tempy = cellCandidData(j,4);
        tempList = SwitchColumn1_2(cat(1,peakGroupStats(cellCandidData(j,1)).PixelList));
        idx =  knnsearch(tempList(:,3),interpZList(round(tempy),tempRowLabel));
        cellCandidData(j,3:5) = tempList(idx,:);
    end
end

% Remove too close candidate pairs (2nd)
for rowNo = 1:3
    row1List = sortrows(cellCandidData(cellCandidData(:,2)==1,3:5),2);
    row2List = sortrows(cellCandidData(cellCandidData(:,2)==2,3:5),2);
    row3List = sortrows(cellCandidData(cellCandidData(:,2)==3,3:5),2);
    allRowList = {row1List;row2List;row3List};
    row1Intervals = diff(row1List);
    row2Intervals = diff(row2List);
    row3Intervals = diff(row3List);
    allRowIntervals = {row1Intervals;row2Intervals;row3Intervals};
    cellCandidCoord = cellCandidData(:,3:5);
    
    for j = 1:size(allRowIntervals{rowNo,1})
        if allRowIntervals{rowNo,1}(j,2) < 5
            xyz1 = allRowList{rowNo,1}(j,:);
            xyz2 = allRowList{rowNo,1}(j+1,:);
            f1 = (cellCandidCoord(:,1) == xyz1(1)) .* (cellCandidCoord(:,2) == xyz1(2)) ...
                .* (cellCandidCoord(:,3) == xyz1(3));
            f2 = (cellCandidCoord(:,1) == xyz2(1)) .* (cellCandidCoord(:,2) == xyz2(2)) ...
                .* (cellCandidCoord(:,3) == xyz2(3));
            if abs(xyz1(1)-xyz2(1))<4
                cellCandidData(f1>0,2) = -1;
                cellCandidData(f2>0,2) = -1;
            else
                cellCandidData(f1>0,2) = 4;
                cellCandidData(f2>0,2) = 4;
            end
        end
    end
end

% Recover removed candidates with high prediction scores (2nd)
f1 = cellCandidData(:,2)==-1;
f2 = cellCandidData(:,6)>0.8;
cellCandidData(f1.*f2>0,2) = 4;

%% Examine gaps within row of estimated outer hair cells 
peakPixelList = GetCoordOfPositivePixels(labeledIm>0);

for rowNo = 1:3
    rowCoordList = sortrows(cellCandidData(cellCandidData(:,2)==rowNo,3:5),2);
    temp = diff(rowCoordList);
    neighborDistList = sqrt(sum(temp(:,1:2).^2,2));
    cellCandidCoord = cellCandidData(:,3:5);
    for j = 1:size(neighborDistList)
        if neighborDistList(j,1) > 11
            leftCellCoord = rowCoordList(j,:);
            rightCellCoord = rowCoordList(j+1,:);            
            tempLeftCellCoord = leftCellCoord;
            neighborDist = norm(leftCellCoord-rightCellCoord); % Estimate the number of lost cell in gap space
            lostCellNo = max(1,round(neighborDist/interpIntervalList(round(leftCellCoord(2)))-1));
            
            if lostCellNo == 1
                estimatedCoord = mean([leftCellCoord;rightCellCoord]);
                predictorIm = ObtainPixelValuesAround2(estimatedCoord,linearizedIm2);
                isGap = classify(net2,predictorIm); % Judge cell or gap
                [~,nearestPeakPoint,nearestDist] = ObtainFeatureQuantities(estimatedCoord,rowNo,selectedPixelList ...
                    ,selectedIndexList,labeledIm,predictionScores2,leftEnd,correlationMatrix); % Search nearest peak point

                if isGap == '-1' % if existence of cell in the gap indicated 
                    if nearestDist < 3.5
                        groupIndex = labeledIm(nearestPeakPoint(1),nearestPeakPoint(2),nearestPeakPoint(3));
                        candidateIndex = selectedIndexList == groupIndex;
                        tempRowNo = cellCandidData(candidateIndex,2);
                        if tempRowNo ~= rowNo
                            cellCandidData(candidateIndex,2) = rowNo;
                            cellCandidData(candidateIndex,3:5) = nearestPeakPoint;
                        end
                    else
                        [peakPixelIndex,nearestDist] = knnsearch(peakPixelList,estimatedCoord); % Search nearest peak point (extended)
                        if nearestDist < 3.5
                            nearestPeakPoint = peakPixelList(peakPixelIndex,:);
                            tempRowNo = labeledIm(nearestPeakPoint(1),nearestPeakPoint(2),nearestPeakPoint(3));
                            cellCandidData = [cellCandidData; double(tempRowNo),rowNo,nearestPeakPoint,NaN];
                        end
                    end
                end
                
            elseif lostCellNo == 2
                vector = (rightCellCoord-leftCellCoord);
                vector = vector/3;
                estimatedCoords = [leftCellCoord + vector;rightCellCoord - vector];
                predictorIm1 = ObtainPixelValuesAround2(estimatedCoords(1,:),linearizedIm2);
                predictorIm2 = ObtainPixelValuesAround2(estimatedCoords(2,:),linearizedIm2);
                isGap = classify(net2,cat(4,predictorIm1, predictorIm2));
                
                for k = 1:2
                    if isGap(k,1) == '-1'
                        [idx, nearestDist] = knnsearch(selectedPixelList,estimatedCoords(k,:));
                        if nearestDist < 3.5
                            nearestPeakPoint = selectedPixelList(idx,:);
                            candidateIndex = selectedIndexList == groupIndex;
                            groupIndex = labeledIm(nearestPeakPoint(1),nearestPeakPoint(2),nearestPeakPoint(3));
                            tempRowNo = cellCandidData(candidateIndex,2);
                            if tempRowNo ~= rowNo
                                cellCandidData(candidateIndex,2) = rowNo;
                                cellCandidData(candidateIndex,3:5) = nearestPeakPoint;
                            end
                        else
                            [peakPixelIndex,nearestDist] = knnsearch(peakPixelList,estimatedCoords(k,:));
                            if nearestDist < 3.5
                                nearestPeakPoint = peakPixelList(peakPixelIndex,:);
                                tempRowNo = labeledIm(nearestPeakPoint(1),nearestPeakPoint(2),nearestPeakPoint(3));
                                cellCandidData = [cellCandidData; double(tempRowNo),rowNo,nearestPeakPoint, NaN];
                            end
                        end
                    end
                end
                
            else % Lost cell number > 2
                for kk = 1:100                    
                    % Estimate cell coordinates in the gap from left to right 
                    row1List = sortrows(cellCandidData(cellCandidData(:,2)==1,3:5),2);
                    row2List = sortrows(cellCandidData(cellCandidData(:,2)==2,3:5),2);
                    row3List = sortrows(cellCandidData(cellCandidData(:,2)==3,3:5),2);                    
                    [idx1] = knnsearch(row1List,tempLeftCellCoord,'k',2);
                    [idx2] = knnsearch(row2List,tempLeftCellCoord,'k',2);
                    [idx3] = knnsearch(row3List,tempLeftCellCoord,'k',2);
                    vector1 = row1List(max(idx1),:)-row1List(min(idx1),:);
                    vector1 = vector1/norm(vector1);
                    vector2 = row2List(max(idx2),:)-row2List(min(idx2),:);
                    vector2 = vector2/norm(vector2);
                    vector3 = row3List(max(idx3),:)-row3List(min(idx3),:);
                    vector3 = vector3/norm(vector3);                    
                    vector4 = rightCellCoord - tempLeftCellCoord;
                    vector4 = vector4/norm(vector4);                    
                    mergedVector = mean([vector1;vector2;vector3])+vector4;
                    mergedVector(3) = 0;
                    mergedVector = mergedVector/norm(mergedVector); % Estimated direction of next cell                    
                    nextCellCoord = tempLeftCellCoord + mergedVector*interpIntervalList(round(leftCellCoord(2)));
                    
                    % Break if distance to right side cell is too small
                    distance1 = norm(rightCellCoord-nextCellCoord);
                    distance2 = rightCellCoord(2)-nextCellCoord(2);
                    if distance2 < interpIntervalList(round(rightCellCoord(2)))/2
                        break
                    elseif distance1 < interpIntervalList(round(rightCellCoord(2)))
                        nextCellCoord = mean([tempLeftCellCoord; rightCellCoord]);
                    end
                    
                    predictorIm = ObtainPixelValuesAround2(nextCellCoord,linearizedIm2);
                    isGap = classify(net2,predictorIm);
                    
                    [~,nearestPeakPoint,nearestDist] = ObtainFeatureQuantities(nextCellCoord,rowNo,selectedPixelList ...
                        ,selectedIndexList,labeledIm,predictionScores2,leftEnd,correlationMatrix); % Search nearest peak point
                    
                    if isGap == '-1'
                        if nearestDist < 3.5
                            groupIndex = labeledIm(nearestPeakPoint(1),nearestPeakPoint(2),nearestPeakPoint(3));
                            candidateIndex = selectedIndexList == groupIndex;
                            tempRowNo = cellCandidData(candidateIndex,2);
                            if tempRowNo ~= rowNo 
                                if tempRowNo == -1 || tempRowNo == 4
                                    cellCandidData(candidateIndex,2) = rowNo;
                                    cellCandidData(candidateIndex,3:5) = nearestPeakPoint;
                                elseif rowNo == 1 
                                    nearestPeakPoint = tempLeftCellCoord;
                                    nearestPeakPoint(1) = nearestPeakPoint(1)-5;
                                elseif rowNo == 2 
                                    cellCandidData(candidateIndex,2) = rowNo;
                                    cellCandidData(candidateIndex,3:5) = nearestPeakPoint;
                                elseif rowNo == 3 
                                    nearestPeakPoint = tempLeftCellCoord; 
                                    nearestPeakPoint(1) = nearestPeakPoint(1)+5;
                                end
                            end
                            nextCellCoord = nearestPeakPoint;
                        else
                            [peakPixelIndex,nearestDist] = knnsearch(peakPixelList,nextCellCoord);
                            if nearestDist < 3.5
                                nearestPeakPoint = peakPixelList(peakPixelIndex,:);
                                tempRowNo = labeledIm(nearestPeakPoint(1),nearestPeakPoint(2),nearestPeakPoint(3));
                                cellCandidData = [cellCandidData; double(tempRowNo),rowNo,nearestPeakPoint,NaN];
                                nextCellCoord = nearestPeakPoint;
                            end                        
                        end
                    end
                    tempLeftCellCoord = nextCellCoord; % Prepare for next search
                end
            end
            
        end
    end
end

%% Check basal end
addCoordList = cell(3,1);
for j = 1:200
    row1List = sortrows([cellCandidData(cellCandidData(:,2)==1,3:5); addCoordList{1,1}],2);
    row2List = sortrows([cellCandidData(cellCandidData(:,2)==2,3:5); addCoordList{2,1}],2);
    row3List = sortrows([cellCandidData(cellCandidData(:,2)==3,3:5); addCoordList{3,1}],2);
    basalEndList = [row1List(end,:);row2List(end,:);row3List(end,:)];
    [~,currentRowNo] = min(basalEndList(:,2));
    tempLeftCellCoord = basalEndList(currentRowNo,:); % Start point for search

    [idx1] = knnsearch(row1List,tempLeftCellCoord,'k',2);
    [idx2] = knnsearch(row2List,tempLeftCellCoord,'k',2);
    [idx3] = knnsearch(row3List,tempLeftCellCoord,'k',2);
    vector1 = row1List(max(idx1),:)-row1List(min(idx1),:);
    vector1 = vector1/norm(vector1);
    vector2 = row2List(max(idx2),:)-row2List(min(idx2),:);
    vector2 = vector2/norm(vector2);
    vector3 = row3List(max(idx3),:)-row3List(min(idx3),:);
    vector3 = vector3/norm(vector3);
    mergedVector = mean([vector1;vector2;vector3;[0 1 0]]);
    mergedVector = mergedVector/norm(mergedVector); % Estimated direction of next cell    
    nextCellCoord = tempLeftCellCoord + mergedVector*7.8;
    
    % Break if the distance to image boundary is too small
    if nextCellCoord(2) > imSize(2)-10
        break
    end
    nearestPeakPoint = round(nextCellCoord);
    pixelValue = linearizedIm2(nearestPeakPoint(1),nearestPeakPoint(2)+5,nearestPeakPoint(3));
    if pixelValue==0
        break
    end
    
    predictorIm = ObtainPixelValuesAround2(nextCellCoord,linearizedIm2);
    isGap = classify(net2,predictorIm);
    
    [~,nearestPeakPoint,nearestDist] = ObtainFeatureQuantities(nextCellCoord,currentRowNo,selectedPixelList,selectedIndexList,labeledIm ...
        ,predictionScores2,leftEnd,correlationMatrix);
    
    if isGap == '-1' && nearestDist < 3 % if there is peak point near the estimated coordinate
        nextCellCoord = nearestPeakPoint;
        groupIndex = labeledIm(nearestPeakPoint(1),nearestPeakPoint(2),nearestPeakPoint(3));
        candidateIndex = selectedIndexList == groupIndex;
        tempRowNo = cellCandidData(candidateIndex,2);   
         
        if tempRowNo == -1 || tempRowNo == 4
            cellCandidData(candidateIndex,2) = currentRowNo;
            cellCandidData(candidateIndex,3:5) = nextCellCoord;
        elseif currentRowNo == 1 % Avoid crossing the other rows
            nextCellCoord = tempLeftCellCoord;
            nextCellCoord(1) = nextCellCoord(1)-5;
            nextCellCoord(2) = nextCellCoord(2)+1;
            addCoordList{currentRowNo,1} = [addCoordList{currentRowNo,1}; nextCellCoord];
        elseif currentRowNo == 2
            cellCandidData(candidateIndex,2) = currentRowNo;
            cellCandidData(candidateIndex,3:5) = nextCellCoord;
        elseif currentRowNo == 3
            cellCandidData(candidateIndex,2) = currentRowNo;
            cellCandidData(candidateIndex,3:5) = nextCellCoord;
        end        
    elseif isGap == '-1' % If there is no peak point nearby, just record the coordinate
        addCoordList{currentRowNo,1} = [addCoordList{currentRowNo,1}; nextCellCoord];
    else
        [peakPixelIndex,nearestDist] = knnsearch(peakPixelList,nextCellCoord);
        if nearestDist < 3 
            nearestPeakPoint = peakPixelList(peakPixelIndex,:);
            tempRowNo = labeledIm(nearestPeakPoint(1),nearestPeakPoint(2),nearestPeakPoint(3));
            cellCandidData = [cellCandidData; double(tempRowNo),currentRowNo,nearestPeakPoint,NaN];
        else
            addCoordList{currentRowNo,1} = [addCoordList{currentRowNo,1}; nextCellCoord];
        end
    end
end

row1List = sortrows(cellCandidData(cellCandidData(:,2)==1,3:5),2);
row2List = sortrows(cellCandidData(cellCandidData(:,2)==2,3:5),2);
row3List = sortrows(cellCandidData(cellCandidData(:,2)==3,3:5),2);
outerHairCells = [row1List; row2List; row3List];

% Add cell candidates labelled with "4" (this label means "pending") for the final check
pendingCoordList = sortrows(cellCandidData(cellCandidData(:,2)==4,3:5),2);
pendingCoordList2 = [];
for k = 1:20 % Remove candidates too close to existing points 
    if isempty(pendingCoordList2)
        [~,nearestDist] = knnsearch(outerHairCells(:,1:2),pendingCoordList(:,1:2));
    else
        [~,nearestDist] = knnsearch([outerHairCells(:,1:2); pendingCoordList2(:,1:2)],pendingCoordList(:,1:2));
    end
    f = (nearestDist>6).*(nearestDist<10);
    pendingCoordList2 = [pendingCoordList2; pendingCoordList(f>0,:)];
end
f = abs(pendingCoordList2(:,3) - interpZList(round(pendingCoordList2(:,2)),3))>5;
pendingCoordList2(f>0,:) = []; % Remove candidates far from average z coordinate

% Add cell candidates at apical end for the final check 
tempList = [];
for k = 1:size(apicalCandidates,1)
    f = cellCandidData(:,1)==apicalCandidates(k,1);
    if cellCandidData(f>0,2)==-1 || cellCandidData(f>0,2)==4
        [~,nearestDist] = knnsearch([outerHairCells(:,1:2); pendingCoordList2(:,1:2)],cellCandidData(f>0,3:4));
        if nearestDist>5
           tempList = [tempList; find(f)];
        else
            nearestDist;
        end
    end
end
apicalCoordList = cellCandidData(tempList,3:5);
f = abs(apicalCoordList(:,3) - mean(interpZList(round(apicalCoordList(:,2)),:),2))>5;
apicalCoordList(f>0,:) = [];
pendingCoordList2 = unique([pendingCoordList2;apicalCoordList],'rows');

% Final check
predictorImList = zeros(69,39,1,size(pendingCoordList2,1));
for i = 1:size(pendingCoordList2,1)
    predictorImList(:,:,:,i) = ObtainPixelValuesAround2(pendingCoordList2(i,:),linearizedIm2);
end
isGap = classify(net2,predictorImList);
pendingCoordList2(isGap=='1',:)=[];

outerHairCells = [sortrows(cellCandidData(cellCandidData(:,2)==1,3:5),2); sortrows(cellCandidData(cellCandidData(:,2) ...
    ==2,3:5),2); sortrows(cellCandidData(cellCandidData(:,2)==3,3:5),2); pendingCoordList2];

markedImOfOHC = DrawMarkedIm(linearizedIm2,round(outerHairCells),1,1);
markedImOfOHC = uint16(markedImOfOHC*1000);
ImWrite3D(markedImOfOHC,[SubFolderPath '\outerHairCells.tif']);
xlswrite([SubFolderPath '\outerHairCells.xlsx'],sortrows(SwitchColumn1_2(round(outerHairCells))),1)
save([SubFolderPath '\data.mat'],'outerHairCells','-append')

disp('outer hair cells estimated!')
