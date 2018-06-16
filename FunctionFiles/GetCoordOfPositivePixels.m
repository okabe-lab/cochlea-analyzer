function xyz = GetCoordOfPositivePixels(binary_image)
%*****************************************************************************
%This function uses "binary_image", binarized image (stack). The output is
%"xyz", coordinates of positive pixels.
%*****************************************************************************
binary_image = binary_image>0;
Isize = size(binary_image);
if size(Isize,2)==2
    [X,Y] = ind2sub([Isize(1),Isize(2)],find(binary_image==1));
    xyz = [X,Y];
    return
end
[X,Y,Z] = ind2sub([Isize(1),Isize(2),Isize(3)],find(binary_image==1));
xyz = [X,Y,Z];
end