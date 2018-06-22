%*****************************************************************************
%This script simulate cell loss of organ of corti with "Neighborhood
%Effect" and "Positional Effect". The script uses MAT-File "analyzeResults.mat", 
%which is created by "main4.m". The script will show three histograms of cluster 
%size and an image indicating similarity between measured and simulated values.
%*****************************************************************************

if not(exist('mainPath','var'))
    mainPath = fileparts(mfilename('fullpath')); % Get fullpath of this script
end
cd(mainPath)

%% Information input
%Specify the full path of a folder containing "analyzeResults.mat".
SubFolderPath = '.\TestData\Result';

%% Load measured sizes of clusters of cell loss
load([SubFolderPath '\analyzeResults.mat'],'totalLossRatio','lossClusterSizes')
Table0 = tabulate(lossClusterSizes);
originalHistogram = Table0(:,1).*Table0(:,2);
originalHistogram(10) = sum(originalHistogram(10:end));
originalHistogram(11:end)=[];
originalHistogram = originalHistogram/sum(originalHistogram);

%% Simulation 1, "Neighborhood neighborEffect" only
Cochleasize = [3 600];
neighborEffect = [0 1 2 3 4 6 10 12];
trialNo = 1000;
errorScores1 = zeros(numel(neighborEffect),1);
% Search value of neighbor effect with minimum error
for i = 1:numel(neighborEffect)
    fprintf('Model 1 [%d / %d] ...\n',i,numel(neighborEffect))
    [totalcellMat] = SimulateNeighborEffect(Cochleasize, totalLossRatio, neighborEffect(i),trialNo);
    Table1 = ComputeSizeDistribution(totalcellMat, totalLossRatio);
    
    model1Histogram = Table1(:,1).*Table1(:,2);
    model1Histogram(10) = sum(model1Histogram(10:end));
    model1Histogram(11:end)=[];
    model1Histogram = model1Histogram/sum(model1Histogram);
    % Error score = RSS / Maximum RSS
    errorScores1(i) = sqrt(sum((model1Histogram - originalHistogram).^2))...
        /sqrt((sum(model1Histogram.^2)-min(model1Histogram)^2+(1-min(model1Histogram))^2));
end

[~,minErrorIdx1] = min(errorScores1); % Recalculation for making graph
[totalcellMat] = SimulateNeighborEffect(Cochleasize, totalLossRatio, neighborEffect(minErrorIdx1),trialNo);
Table1 = ComputeSizeDistribution(totalcellMat, totalLossRatio);
model1Histogram = Table1(:,1).*Table1(:,2);
model1Histogram(10) = sum(model1Histogram(10:end));
model1Histogram(11:end)=[];
model1Histogram = model1Histogram/sum(model1Histogram);

figure
bar([originalHistogram(1:10) model1Histogram(1:10)])
ax = gca;
ax.FontSize = 14;
legend({'Measured','Simulation 1'})
xlabel('Number in cluster')
ylabel('Probability')
title(sprintf('Ratio= %.3f, Neighbor= %d, ErrScore= %.2f',totalLossRatio,neighborEffect(minErrorIdx1),min(errorScores1(:))))
xticklabels({'1','2','3','4','5','6','7','8','9','ÅÜ10'})

%% Simulation 2, "Positional neighborEffect" only
powerInd = [0.01 1 3 5 10];
gaussInd = [1 2 4 6];
errorScores2 = zeros(numel(powerInd),numel(gaussInd));
count = 0;
totalCount = numel(powerInd)*numel(gaussInd);
for i = 1:numel(powerInd)
    for j = 1:numel(gaussInd)
        count = count+1;
        fprintf('Model 2 [%d / %d] ...\n',count,totalCount);
        [totalcellMat] = SimulatePositionalEffect(Cochleasize, totalLossRatio, powerInd(i) ...
            ,gaussInd(j),trialNo);
        Table2 = ComputeSizeDistribution(totalcellMat, totalLossRatio);
        
        model2Histogram = Table2(:,1).*Table2(:,2);
        model2Histogram(10) = sum(model2Histogram(10:end));
        model2Histogram(11:end)=[];
        model2Histogram = model2Histogram/sum(model2Histogram);
        % Error score = RSS / Maximum RSS
        errorScores2(i,j) = sqrt(sum((model2Histogram - originalHistogram).^2))...
            /sqrt((sum(model2Histogram.^2)-min(model2Histogram)^2+(1-min(model2Histogram))^2));
    end
end

[~,minErrorIdx2] = min(errorScores2(:));
[i,j] = ind2sub(size(errorScores2),minErrorIdx2);
[totalcellMat] = SimulatePositionalEffect(Cochleasize, totalLossRatio, powerInd(i),gaussInd(j),trialNo);
Table2 = ComputeSizeDistribution(totalcellMat, totalLossRatio);

model2Histogram = Table2(:,1).*Table2(:,2);
model2Histogram(10) = sum(model2Histogram(10:end));
model2Histogram(11:end)=[];
model2Histogram = model2Histogram/sum(model2Histogram);

figure
bar([originalHistogram(1:10) model2Histogram(1:10)])
ax = gca;
ax.FontSize = 14;
legend({'Measured','Simulation 2'})
xlabel('Number in cluster')
ylabel('Probability')
title(sprintf('Ratio= %.3f, Power= %.3f, Gauss= %d, ErrScore= %.2f',totalLossRatio,powerInd(i) ...
    ,gaussInd(j),min(errorScores2(:))))
xticklabels({'1','2','3','4','5','6','7','8','9','ÅÜ10'})

%% Simulation 3, combination of two effects
trialNo = 500;
Cochleasize = [3 600];
neighborEffect = [0 2 4 8];
positionEffect = [0 4.5 6 8];
errorScores3 = zeros(numel(neighborEffect),numel(positionEffect));
count = 0;
totalCount = numel(neighborEffect)*numel(positionEffect);
for i = 1:numel(neighborEffect)
    for j = 1:numel(positionEffect)
        count = count + 1;
        fprintf('Model 1+2 [%d / %d] ...\n',count,totalCount);
        [totalcellMat] = SimulateCombinationEffect(Cochleasize, totalLossRatio, neighborEffect(i),positionEffect(j) ...
            ,trialNo);
        Table4 = ComputeSizeDistribution(permute(totalcellMat,[2 1 3]), totalLossRatio);
        
        model3Histogram = Table4(:,1).*Table4(:,2);
        model3Histogram(10) = sum(model3Histogram(10:end));
        model3Histogram(11:end)=[];
        model3Histogram = model3Histogram/sum(model3Histogram);
        % Error score = RSS / Maximum RSS
        errorScores3(i,j) = sqrt(sum((model3Histogram - originalHistogram).^2))...
            /sqrt((sum(model3Histogram.^2)-min(model3Histogram)^2+(1-min(model3Histogram))^2));
    end
end

[~,minErrorIdx3] = min(errorScores3(:));
[i,j] = ind2sub(size(errorScores3),minErrorIdx3);
[totalcellMat] = SimulateCombinationEffect(Cochleasize, totalLossRatio, neighborEffect(i), positionEffect(j),trialNo);
Table4 = ComputeSizeDistribution(totalcellMat, totalLossRatio);

model3Histogram = Table4(:,1).*Table4(:,2);
model3Histogram(10) = sum(model3Histogram(10:end));
model3Histogram(11:end)=[];
model3Histogram = model3Histogram/sum(model3Histogram);

% Weighted average of top 3 combinations
[~,idx] = sort(errorScores3(:));
temp = -errorScores3-(-errorScores3(idx(end-3)));
temp(idx(4:end))=0;
temp = temp.^2;
tempx = temp.*repmat([1;2;3;4],[1,4]);
centerX = sum(tempx(:))/sum(temp(:));
tempy = temp.*repmat([1 2 3 4],[4,1]);
centerY = sum(tempy(:))/sum(temp(:));

[~,minErrorIdx3] = min(errorScores3(:));
[i,j] = ind2sub(size(errorScores3),minErrorIdx3);

figure
bar([originalHistogram(1:10) model3Histogram(1:10)])
ax = gca;
ax.FontSize = 14;
legend({'Measured','Simulation 3'})
xlabel('Number in cluster')
ylabel('Probability')
title(sprintf('Ratio= %.3f, Neighbor = %d, Area = %d, ErrScore= %.2f',totalLossRatio,neighborEffect(i) ...
    ,positionEffect(j),min(errorScores3(:))))
xticklabels({'1','2','3','4','5','6','7','8','9','ÅÜ10'})

figure
image(flipud(-errorScores3),'CDataMapping','scaled')
hold on
colorbar
xticks(1:4)
yticks(1:4)
xticklabels(0:3)
yticklabels(flip(0:3))
ax = gca;
ax.YDir = 'reverse';
scatter(centerY,4-centerX+1,'filled','r')
hold off
xlabel('Positional')
ylabel('Neighborhood')
ax.FontSize = 14;

disp('main6.m completed!')