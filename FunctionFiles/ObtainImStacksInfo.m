function [TiffFileNames, xyStageCoord_micron] = ObtainImStacksInfo(folderPath ...
    ,excelFileName)
%*****************************************************************************
%This function uses "folderPath", path of folder for search image stacks
%and excel file, and "excelFileName", excel file containing imaging
%positions on XY scanning stage. The outputs are "fileNameList", names of 
%image stacks and "xyStageCoord_micron", xy(z) coordinates of xy scanning 
%stage.
%*****************************************************************************

% Get file names of image stacks and read excel file with
% coordinates on XY scanning stage for each image stacks.
folderInfo = dir(folderPath);
fileNames = cat(1,{folderInfo.name})';
isTiffFile = FindStringPattern(fileNames,'.tif');
TiffFileNames = fileNames(isTiffFile > 0);

temp = xlsread([folderPath '\' excelFileName]);
xyStageCoord_micron = temp(:,2:3);