%*****************************************************************************
%This script make heatmap image indicating degrees of cell loss. The script
%uses multiple MAT-file "analyzeResults.mat", which are created by the script 
%"main4.m". 
%*****************************************************************************

if not(exist('mainPath','var'))
    mainPath = fileparts(mfilename('fullpath')); % Get fullpath of this script
end
cd(mainPath)

%% Information input
% Specify cell arrays of paths of ImFolderPathList containing "analyzeResults.mat".
% For example, {'G:\Cochlear\No1\Result','G:\Cochlear\No2\Result'}
ImFolderPathList = {'.\TestData\Result','.\TestData\Result','.\TestData\Result','.\TestData\Result'}; 

%% Draw heatmap image of degree of cell loss
void_values50s = cell(numel(ImFolderPathList),1);
ratios = zeros(numel(ImFolderPathList),1);

for j = 1:numel(ImFolderPathList)
    load([ImFolderPathList{j} '\analyzeResults.mat'],'lostCellNoList50','totalLossRatio');
    void_values50s{j} = lostCellNoList50;
    ratios(j,1) = totalLossRatio;
end

void_values1 = cat(2,void_values50s{:});

[~,ind] = sort(ratios,'descend');
values1 = void_values1(:,ind);
values3 = ratios(ind);

values = values1';

figure
xvalues = num2cell(0:50:5950);
yvalues = [];
h = heatmap(values,'Colormap',parula,'ColorLimits',[0 6]);
ax = gca;
ax.FontSize = 14;

disp('main5.m completed!')