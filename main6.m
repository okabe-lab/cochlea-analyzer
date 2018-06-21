%*****************************************************************************
%This script simulate cell loss of organ of corti with "Neighborhood
%Effect" and "Positional Effect". The script uses MAT-File "result2.mat",
%which is created by "main4.m". The script will show three histograms of
%cluster size and an image indicating similarity between measured and
%simulated values.
%*****************************************************************************

if not(exist('mainPath','var'))
    mainPath = fileparts(mfilename('fullpath')); % Get fullpath of this script
end
cd(mainPath)

%% Information input
%Specify the full path of a folder containing "analyzeResults.mat".
SubFolderPath = '.\TestData\Result';

%% Load measured sizes of clusters of cell loss
load([SubFolderPath '\analyzeResults.mat'],'Ratio2','enum2')
enum = enum2;
Ratio = Ratio2;
Table0 = tabulate(enum);
T0 = Table0(:,1).*Table0(:,2);
T0(10) = sum(T0(10:end));
T0(11:end)=[];
T0 = T0/sum(T0);

%% Simulation 1, "Neighborhood Effect" only
Cochleasize = [3 600];
effect = [0 1 2 3 4 6 10 12];
trialNo = 1000;
er1_2 = zeros(numel(effect),1);
for i = 1:numel(effect)
    fprintf('Model 1 [%d / %d] ...\n',i,numel(effect))
    [totalcellMat] = SimulateNeighborEffect(Cochleasize, Ratio, effect(i),trialNo);
    Table1 = ComputeSizeDistribution(totalcellMat, Ratio);
    
    T1 = Table1(:,1).*Table1(:,2);
    T1(10) = sum(T1(10:end));
    T1(11:end)=[];
    T1 = T1/sum(T1);
    er1_2(i) = sqrt(sum((T1 - T0).^2))/sqrt((sum(T1.^2)-min(T1)^2+(1-min(T1))^2));
end

[~,idx1_2] = min(er1_2);
[totalcellMat] = SimulateNeighborEffect(Cochleasize, Ratio, effect(idx1_2),trialNo);
Table1 = ComputeSizeDistribution(totalcellMat, Ratio);
T1 = Table1(:,1).*Table1(:,2);
T1(10) = sum(T1(10:end));
T1(11:end)=[];
T1 = T1/sum(T1);

ef2 = effect(idx1_2);

figure
bar([T0(1:10) T1(1:10)])
ax = gca;
ax.FontSize = 14;
legend({'Measured','Simulation 1'})
xlabel('Number in cluster')
ylabel('Probability')
title(sprintf('Ratio= %.3f, Neighbor= %d, ErrScore= %.1f',Ratio,effect(idx1_2),min(er1_2(:))*100))
xticklabels({'1','2','3','4','5','6','7','8','9','ÅÜ10'})

%% Simulation 2, "Positional Effect" only
powerInd = [0.01 1 3 5 10];
gaussInd = [1 2 4 6];
er2_2 = zeros(numel(powerInd),numel(gaussInd));
count = 0;
totalCount = numel(powerInd)*numel(gaussInd);
for i = 1:numel(powerInd)
    for j = 1:numel(gaussInd)
        count = count+1;
        fprintf('Model 2 [%d / %d] ...\n',count,totalCount);
        [totalcellMat] = SimulatePositionalEffect(Cochleasize, Ratio, powerInd(i) ...
            ,gaussInd(j),trialNo);
        Table2 = ComputeSizeDistribution(totalcellMat, Ratio);
        
        T2 = Table2(:,1).*Table2(:,2);
        T2(10) = sum(T2(10:end));
        T2(11:end)=[];
        T2 = T2/sum(T2);
        er2_2(i,j) = sqrt(sum((T2 - T0).^2))/sqrt((sum(T2.^2)-min(T2)^2+(1-min(T2))^2));
    end
end

[~,idx2_2] = min(er2_2(:));
[i,j] = ind2sub(size(er2_2),idx2_2);
[totalcellMat] = SimulatePositionalEffect(Cochleasize, Ratio, powerInd(i),gaussInd(j),trialNo);
Table2 = ComputeSizeDistribution(totalcellMat, Ratio);

T2 = Table2(:,1).*Table2(:,2);
T2(10) = sum(T2(10:end));
T2(11:end)=[];
T2 = T2/sum(T2);

pow2 = powerInd(i);
gau2 = gaussInd(j);

figure
bar([T0(1:10) T2(1:10)])
ax = gca;
ax.FontSize = 14;
legend({'Measured','Simulation 2'})
xlabel('Number in cluster')
ylabel('Probability')
title(sprintf('Ratio= %.3f, Power= %.3f, Gauss= %d, ErrScore= %.1f',Ratio,powerInd(i) ...
    ,gaussInd(j),min(er2_2(:))*100))
xticklabels({'1','2','3','4','5','6','7','8','9','ÅÜ10'})

%% Simulation 3, combination of two effects
trialNo = 500;
Cochleasize = [3 600];
effect = [0 2 4 8];
areaInd = [0 4.5 6 8];
er4_2 = zeros(numel(effect),numel(areaInd));
count = 0;
totalCount = numel(effect)*numel(areaInd);
for i = 1:numel(effect)
    for j = 1:numel(areaInd)
        count = count + 1;
        fprintf('Model 1+2 [%d / %d] ...\n',count,totalCount);
        [totalcellMat] = SimulateCombinationEffect(Cochleasize, Ratio, effect(i),areaInd(j) ...
            ,trialNo);
        Table4 = ComputeSizeDistribution(permute(totalcellMat,[2 1 3]), Ratio);
        
        T4 = Table4(:,1).*Table4(:,2);
        T4(10) = sum(T4(10:end));
        T4(11:end)=[];
        T4 = T4/sum(T4);
        er4_2(i,j) = sqrt(sum((T4 - T0).^2))/sqrt((sum(T4.^2)-min(T4)^2+(1-min(T4))^2));
    end
end

[~,idx2_2] = min(er4_2(:));
[i,j] = ind2sub(size(er4_2),idx2_2);
[totalcellMat] = SimulateCombinationEffect(Cochleasize, Ratio, effect(i), areaInd(j),trialNo);
Table4 = ComputeSizeDistribution(totalcellMat, Ratio);

T4 = Table4(:,1).*Table4(:,2);
T4(10) = sum(T4(10:end));
T4(11:end)=[];
T4 = T4/sum(T4);

[~,idx] = sort(er4_2(:));
temp = -er4_2-(-er4_2(idx(end-3)));
temp(idx(4:end))=0;
temp = temp.^2;
tempx = temp.*repmat([1;2;3;4],[1,4]);
centerx = sum(tempx(:))/sum(temp(:));
tempy = temp.*repmat([1 2 3 4],[4,1]);
centery = sum(tempy(:))/sum(temp(:));

[~,idx2_2] = min(er4_2(:));
[i,j] = ind2sub(size(er4_2),idx2_2);

figure
bar([T0(1:10) T4(1:10)])
ax = gca;
ax.FontSize = 14;
legend({'Measured','Simulation 3'})
xlabel('Number in cluster')
ylabel('Probability')
title(sprintf('Ratio= %.3f, Neighbor = %d, Area = %d, ErrScore= %.1f',Ratio,effect(i) ...
    ,areaInd(j),min(er4_2(:))*100))
xticklabels({'1','2','3','4','5','6','7','8','9','ÅÜ10'})

figure
image(flipud(-er4_2*100),'CDataMapping','scaled')
hold on
colorbar
xticks(1:4)
yticks(1:4)
xticklabels(0:3)
yticklabels(flip(0:3))
ax = gca;
ax.YDir = 'reverse';
scatter(centery,4-centerx+1,'filled','r')
hold off
xlabel('Positional')
ylabel('Neighborhood')
ax.FontSize = 14;

disp('main6.m completed!')