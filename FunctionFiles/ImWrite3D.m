function ImWrite3D(image,name,varargin)
%*****************************************************************************
%This function uses "image", image stack to export, "name", file name. 
%*****************************************************************************

%% Export image stack
if nargin == 3
    name = [cell2mat(varargin) '\' name];
end

list = dir(name);
if not(isempty(list))
    delete(name);
end
num_images = size(image,3);
    for i = 1:num_images
        imwrite(image(:,:,i), name, 'WriteMode', 'append', 'Compression', 'none');
    end