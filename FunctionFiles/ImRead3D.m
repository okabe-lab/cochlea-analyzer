function imStack = ImRead3D(imStackName,folderPath)
%*****************************************************************************
%This function uses "imStackName", name of image stack and "folderPath",
%path of folder containing image stack. The output is "imStack", 3D image
%data.
%*****************************************************************************

% Read image stack
if nargin == 1
    info = imfinfo(imStackName); num_images2 = numel(info);
    Isize2 = size(imread(imStackName));
    imStack = zeros(Isize2(1),Isize2(2),num_images2);
    for i = 1:num_images2
        imStack(:,:,i) = imread(imStackName, i, 'Info', info);
    end
    
else
    path = [folderPath '\' imStackName];
    info = imfinfo(path); num_images2 = numel(info);
    Isize2 = size(imread(path));
    imStack = zeros(Isize2(1),Isize2(2),num_images2);
    for i = 1:num_images2
        imStack(:,:,i) = imread(path, i, 'Info', info);
    end
end