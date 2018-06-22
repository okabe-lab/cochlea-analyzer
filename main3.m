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
[predictor1,~,predictionScores1,S, L, outC] = PredictOuterHairCells(linearizedIm2,innerCells2,Mdl);
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
    temp = SwitchColumn1_2(S(selectedIndexList(j)).PixelList);
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
selectedPixelList = SwitchColumn1_2(cat(1,S(selectedIndexList).PixelList));
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
            tempCoord = round(tempData(3:5));
            if tempCoord(1)+5>imSize(1)
                cellCandidData(j,2) = -1;
            else              
                pixelValue1 = linearizedIm2(tempCoord(1),tempCoord(2)+10,tempCoord(3));
                pixelValue2 = linearizedIm2(tempCoord(1)-5,tempCoord(2),tempCoord(3));
                pixelValue3 = linearizedIm2(tempCoord(1)+5,tempCoord(2),tempCoord(3));
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
                candidPixelList = SwitchColumn1_2(cat(1,S(idxList).PixelList));
                [idxs3,d] = knnsearch(candidPixelList,gapCoords);
                candidPixelList2 = candidPixelList(idxs3(d<3,:),:);
                if ~isempty(candidPixelList2)
                    linearIndice = sub2ind(imSize,candidPixelList2(:,1),candidPixelList2(:,2),candidPixelList2(:,3));
                    groupNo1 = L(linearIndice);
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
        tempList = SwitchColumn1_2(cat(1,S(cellCandidData(j,1)).PixelList));
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
o_del = cell(3,1);
delcand = cell(3,1);
addcand = cell(3,1);
addlabel = [];
xyzlistL = GetCoordOfPositivePixels(L>0);

for rowNo = 1:3
    rowCoordList = sortrows(cellCandidData(cellCandidData(:,2)==rowNo,3:5),2);
    temp = diff(rowCoordList);
    difnorm = sqrt(sum(temp(:,1:2).^2,2));
    cellCandidCoord = cellCandidData(:,3:5);
    for j = 1:size(difnorm)
        if difnorm(j,1) > 11
            xyz1 = rowCoordList(j,:);
            xyz2 = rowCoordList(j+1,:);
            
            xyz1_3 = xyz1;
            xyz2_3 = xyz2;
            mindist = norm(xyz1-xyz2);
            
            xyz1_4 = xyz1_3;
            xyz2_4 = xyz2_3;
            dnum = max(1,round(mindist/interpIntervalList(round(xyz1(2)))-1));
            
            if dnum == 1
                xyz3 = mean([xyz1_4;xyz2_4]);
                predictorIm = ObtainPixelValuesAround2(xyz3,linearizedIm2);
                temppred = classify(net2,predictorIm);
                delcand{rowNo,1} = [delcand{rowNo,1}; xyz3];
                [param2,tempCoord,d] = ObtainFeatureQuantities(xyz3,rowNo,selectedPixelList ...
                    ,selectedIndexList,L,predictionScores2,leftEnd,outC);

                if temppred == '1'
                    o_del{rowNo,1} = [o_del{rowNo,1}; xyz3];
                else
                    if d < 3.5
                        tempIdx = L(tempCoord(1),tempCoord(2),tempCoord(3));
                        templabel = cellCandidData(selectedIndexList == tempIdx,2);
                        tempidx2 = find(selectedIndexList == tempIdx);
                        if templabel ~= rowNo
                            cellCandidData(tempidx2,2) = rowNo;
                            cellCandidData(tempidx2,3:5) = tempCoord;
                        end
                    else
                        [Lidx,Ld] = knnsearch(xyzlistL,xyz3);
                        if Ld < 3.5
                            tempCoord = xyzlistL(Lidx,:);
                            templabel = L(tempCoord(1),tempCoord(2),tempCoord(3));
                            cellCandidData = [cellCandidData; double(templabel),rowNo,tempCoord,NaN];
                            addlabel = [addlabel; [double(templabel) xyz3 predictionScores1(templabel,2)]];
                        end
                        addcand{rowNo,1} = [addcand{rowNo,1}; xyz3];
                    end
                end
            elseif dnum == 2
                vector = (xyz2_4-xyz1_4);
                vector = vector/3;
                xyz3 = [xyz1_4 + vector;xyz2_4 - vector];
                param1_1 = ObtainPixelValuesAround2(xyz3(1,:),linearizedIm2);
                param1_2 = ObtainPixelValuesAround2(xyz3(2,:),linearizedIm2);
                temppred = classify(net2,cat(4,param1_1, param1_2));
                
                delcand{rowNo,1} = [delcand{rowNo,1}; xyz3];
                
                param2_1 = ObtainFeatureQuantities(xyz3(1,:),rowNo,selectedPixelList,selectedIndexList,L ...
                    ,predictionScores2,leftEnd,outC);
                param2_2 = ObtainFeatureQuantities(xyz3(2,:),rowNo,selectedPixelList,selectedIndexList,L ...
                    ,predictionScores2,leftEnd,outC);
                
                for k = 1:2
                    if temppred(k,1) == '1'
                        o_del{rowNo,1} = [o_del{rowNo,1}; xyz3(k,:)];
                    else
                        [idx, d] = knnsearch(selectedPixelList,xyz3(k,:));
                        if d < 3.5
                            tempCoord = selectedPixelList(idx,:);
                            tempIdx = L(tempCoord(1),tempCoord(2),tempCoord(3));
                            tempidx2 = find(selectedIndexList == tempIdx);
                            templabel = cellCandidData(tempidx2,2);
                            if templabel ~= rowNo
                                cellCandidData(tempidx2,2) = rowNo;
                                cellCandidData(tempidx2,3:5) = tempCoord;
                            end
                        else
                            [Lidx,Ld] = knnsearch(xyzlistL,xyz3(k,:));
                            if Ld < 3.5
                                tempCoord = xyzlistL(Lidx,:);
                                templabel = L(tempCoord(1),tempCoord(2),tempCoord(3));
                                cellCandidData = [cellCandidData; double(templabel),rowNo,tempCoord, NaN];
                                addlabel = [addlabel; [double(templabel) xyz3(k,:) ...
                                    predictionScores1(templabel,2)]];
                            end
                            addcand{rowNo,1} = [addcand{rowNo,1}; xyz3(k,:)];
                        end
                    end
                end
            else 
                for kk = 1:100
                    lastxyz = xyz1_3;

                    row1List = sortrows(cellCandidData(cellCandidData(:,2)==1,3:5),2);
                    row2List = sortrows(cellCandidData(cellCandidData(:,2)==2,3:5),2);
                    row3List = sortrows(cellCandidData(cellCandidData(:,2)==3,3:5),2);
                    
                    [idx1] = knnsearch(row1List,lastxyz,'k',2);
                    [idx2] = knnsearch(row2List,lastxyz,'k',2);
                    [idx3] = knnsearch(row3List,lastxyz,'k',2);
                    vector1 = row1List(max(idx1),:)-row1List(min(idx1),:);
                    vector1 = vector1/norm(vector1);
                    vector2 = row2List(max(idx2),:)-row2List(min(idx2),:);
                    vector2 = vector2/norm(vector2);
                    vector3 = row3List(max(idx3),:)-row3List(min(idx3),:);
                    vector3 = vector3/norm(vector3);
                    
                    vector4 = xyz2_3 - lastxyz;
                    vector4 = vector4/norm(vector4);
                    
                    vector = mean([vector1;vector2;vector3])+vector4;
                    vector(3) = 0;
                    vector = vector/norm(vector);
                    
                    nextxyz = lastxyz + vector*interpIntervalList(round(xyz1(2)));
                    d1 = norm(xyz2_3-nextxyz);
                    d2 = xyz2_3(2)-nextxyz(2);
                    if d2 < interpIntervalList(round(xyz2(2)))/2
                        break
                    elseif d1 < interpIntervalList(round(xyz2(2)))
                        nextxyz = mean([lastxyz; xyz2_3]);
                    end
                    
                    rparam2 = ObtainPixelValuesAround2(nextxyz,linearizedIm2);
                    temppred = classify(net2,rparam2);
                    
                    [rparam,tempCoord,d] = ObtainFeatureQuantities(nextxyz,rowNo,selectedPixelList ...
                        ,selectedIndexList,L,predictionScores2,leftEnd,outC);
                    
                    delcand{rowNo,1} = [delcand{rowNo,1}; nextxyz];
                    if temppred == '1'
                        o_del{rowNo,1} = [o_del{rowNo,1}; nextxyz];
                    else
                        if d < 3.5
                            tempIdx = L(tempCoord(1),tempCoord(2),tempCoord(3));
                            templabel = cellCandidData(selectedIndexList == tempIdx,2);
                            tempidx2 = find(selectedIndexList == tempIdx);
                            if templabel ~= rowNo 
                                if templabel == -1 || templabel == 4
                                    cellCandidData(tempidx2,2) = rowNo;
                                    cellCandidData(tempidx2,3:5) = tempCoord;
                                elseif rowNo == 1 
                                    tempCoord = lastxyz;
                                    tempCoord(1) = tempCoord(1)-5;
                                elseif rowNo == 2 
                                    cellCandidData(tempidx2,2) = rowNo;
                                    cellCandidData(tempidx2,3:5) = tempCoord;
                                elseif rowNo == 3 
                                    tempCoord = lastxyz; 
                                    tempCoord(1) = tempCoord(1)+5;
                                end
                            end
                            nextxyz = tempCoord;
                        else
                            [Lidx,Ld] = knnsearch(xyzlistL,nextxyz);
                            if Ld < 3.5
                                tempCoord = xyzlistL(Lidx,:);
                                templabel = L(tempCoord(1),tempCoord(2),tempCoord(3));
                                cellCandidData = [cellCandidData; double(templabel),rowNo,tempCoord,NaN];
                                addlabel = [addlabel; [double(templabel) nextxyz ...
                                    predictionScores1(templabel,2)]];
                                nextxyz = tempCoord;
                            end
                            addcand{rowNo,1} = [addcand{rowNo,1}; nextxyz];                          
                        end
                    end
                    xyz1_3 = nextxyz;
                end
            end
            
        end
    end
end

%% Examine the end ot rows of estimated outer hair cells 
tempadd = cell(3,1);
for j = 1:200
    row1List = sortrows([cellCandidData(cellCandidData(:,2)==1,3:5); tempadd{1,1}],2);
    row2List = sortrows([cellCandidData(cellCandidData(:,2)==2,3:5); tempadd{2,1}],2);
    row3List = sortrows([cellCandidData(cellCandidData(:,2)==3,3:5); tempadd{3,1}],2);
    out_lasts = [row1List(end,:);row2List(end,:);row3List(end,:)];
    [~,lineidx] = min(out_lasts(:,2));
    lastxyz = out_lasts(lineidx,:);

    [idx1] = knnsearch(row1List,lastxyz,'k',2);
    [idx2] = knnsearch(row2List,lastxyz,'k',2);
    [idx3] = knnsearch(row3List,lastxyz,'k',2);
    vector1 = row1List(max(idx1),:)-row1List(min(idx1),:);
    vector1 = vector1/norm(vector1);
    vector2 = row2List(max(idx2),:)-row2List(min(idx2),:);
    vector2 = vector2/norm(vector2);
    vector3 = row3List(max(idx3),:)-row3List(min(idx3),:);
    vector3 = vector3/norm(vector3);
    vector = mean([vector1;vector2;vector3;[0 1 0]]);
    vector = vector/norm(vector);
    
    nextxyz = lastxyz + vector*7.8;

    if nextxyz(2) > imSize(2)-10
        break
    end
    tempCoord = round(nextxyz);
    tempint = linearizedIm2(tempCoord(1),tempCoord(2)+5,tempCoord(3));
    if tempint==0
        break
    end
    
    rparam2 = ObtainPixelValuesAround2(nextxyz,linearizedIm2);
    temppred = classify(net2,rparam2);
    
    [rparam,xyz4,d] = ObtainFeatureQuantities(nextxyz,lineidx,selectedPixelList,selectedIndexList,L ...
        ,predictionScores2,leftEnd,outC);
    
    delcand{lineidx,1} = [delcand{lineidx,1}; nextxyz];
    
    if temppred ~= '1' && d < 3
        nextxyz = xyz4;

        tempIdx = L(xyz4(1),xyz4(2),xyz4(3));
        templabel = cellCandidData(selectedIndexList == tempIdx,2);
        tempidx2 = find(selectedIndexList == tempIdx);
        
        if templabel ~= lineidx
            if templabel == -1 || templabel == 4
                cellCandidData(tempidx2,2) = lineidx;
                cellCandidData(tempidx2,3:5) = nextxyz;
            elseif lineidx == 1 
                nextxyz = lastxyz; 
                nextxyz(1) = nextxyz(1)-5;
                nextxyz(2) = nextxyz(2)+1;
                tempadd{lineidx,1} = [tempadd{lineidx,1}; nextxyz];
            elseif lineidx == 2 
                cellCandidData(tempidx2,2) = lineidx;
                cellCandidData(tempidx2,3:5) = nextxyz;
            elseif lineidx == 3 
                cellCandidData(tempidx2,2) = lineidx;
                cellCandidData(tempidx2,3:5) = nextxyz;
            end
        end
        
    elseif temppred == '1' 
        o_del{rowNo,1} = [o_del{rowNo,1}; nextxyz];
        tempadd{lineidx,1} = [tempadd{lineidx,1}; nextxyz];
    else
        [Lidx,Ld] = knnsearch(xyzlistL,nextxyz);
        if Ld < 3 
            tempCoord = xyzlistL(Lidx,:);
            templabel = L(tempCoord(1),tempCoord(2),tempCoord(3));
            cellCandidData = [cellCandidData; double(templabel),lineidx,tempCoord,NaN];
            addlabel = [addlabel; [double(templabel) tempCoord predictionScores1(templabel,2)]];
        else
            tempadd{lineidx,1} = [tempadd{lineidx,1}; nextxyz];
        end
    end
end

row1List = sortrows(cellCandidData(cellCandidData(:,2)==1,3:5),2);
row2List = sortrows(cellCandidData(cellCandidData(:,2)==2,3:5),2);
row3List = sortrows(cellCandidData(cellCandidData(:,2)==3,3:5),2);
outerHairCells = [row1List; row2List; row3List];
out4 = sortrows(cellCandidData(cellCandidData(:,2)==4,3:5),2);
delcandlist = [cat(1,delcand{:}); out4];

add4 = [];
for k = 1:20
    if isempty(add4)
        [~,d] = knnsearch(outerHairCells(:,1:2),out4(:,1:2));
    else
        [~,d] = knnsearch([outerHairCells(:,1:2); add4(:,1:2)],out4(:,1:2));
    end
    f = (d>6).*(d<10);
    add4 = [add4; out4(f>0,:)];
end
f = abs(add4(:,3) - interpZList(round(add4(:,2)),3))>5;
add4(f>0,:) = [];

temp = [];
for k = 1:size(apicalCandidates,1)
    f = cellCandidData(:,1)==apicalCandidates(k,1);
    if cellCandidData(f>0,2)==-1 || cellCandidData(f>0,2)==4
        [~,d] = knnsearch([outerHairCells(:,1:2); add4(:,1:2)],cellCandidData(find(f),3:4));
        if d>5
           temp = [temp; find(f)];
        else
            d;
        end
    end
end

add5 = cellCandidData(temp,3:5);
f = abs(add5(:,3) - mean(interpZList(round(add5(:,2)),:),2))>5;
add5(f>0,:) = [];
add4 = unique([add4;add5],'rows');

tempparam = zeros(69,39,1,size(add4,1));
for i = 1:size(add4,1)
    tempparam(:,:,:,i) = ObtainPixelValuesAround2(add4(i,:),linearizedIm2);
end
temppred = classify(net2,tempparam);
add4(temppred=='1',:)=[];

outlist = [sortrows(cellCandidData(cellCandidData(:,2)==1,3:5),2); sortrows(cellCandidData(cellCandidData(:,2) ...
    ==2,3:5),2); sortrows(cellCandidData(cellCandidData(:,2)==3,3:5),2); add4];

incim = DrawMarkedIm(linearizedIm2,round(outlist),1,1);
incim = uint16(incim*1000);
ImWrite3D(incim,[SubFolderPath '\outerHairCells.tif']);
xlswrite([SubFolderPath '\outerHairCells.xlsx'],sortrows(SwitchColumn1_2(round(outlist))),1)
save([SubFolderPath '\data.mat'],'outlist','-append')
disp('outer hair cells estimated!')
