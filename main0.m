%*****************************************************************************
%This script display version information of your Matlab and unzip test data.
%*****************************************************************************

mainPath = fileparts(mfilename('fullpath')); % Get fullpath of this script
cd(mainPath)
addpath('.\FunctionFiles') 

% Display version information
disp('>>>'); disp('>>>')
disp('This program needs:')
disp('   - Matlab (R2017b or newer version),')
disp('   - Image Processing Toolbox,')
disp('   - Statistics and Machine Learning Toolbox,')
disp('   - Neural Network Toolbox.')
disp('Please check version information of your Matlab.')
disp('>>>'); disp('>>>')
ver
disp('>>>'); disp('>>>')

% Unzip test image data
folderInfo = dir('.\TestData');
fileNames = cat(1,{folderInfo.name})';
isTiffFile = FindStringPattern(fileNames,'.tif');
TiffFileNames = fileNames(isTiffFile > 0);

if size(TiffFileNames,1)~=16
    disp('Unzipping image files...')
    for i = 1:16
        zipFileName = sprintf('416_0%02d.zip',i);
        unzip(['.\TestData\' zipFileName],'.\TestData')
    end
end
disp('main0.m completed!')