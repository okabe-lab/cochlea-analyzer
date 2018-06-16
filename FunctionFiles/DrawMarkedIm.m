function newim = DrawMarkedIm(image,vertex,w,d)
%*****************************************************************************
%This function uses "image", reference image stack to get size information,
%"vertex", xyz coordinates, and "w" and "d", width and depth of marks, 
%respectively. The output is "newim", image with marks. 
%*****************************************************************************

%% Create an image stack with marks indicating specified coordinates
    XYZ = vertex;
    Isize = size(image);
    binary_image = zeros(Isize(1),Isize(2),Isize(3));
    w = round(w/2);

    count = 0;
    for i = 1:size(XYZ,1)
        count = count+1;
        xx = round(XYZ(i,1)); yy = round(XYZ(i,2)); zz = round(XYZ(i,3));
        for Z = zz-d:zz+d
            for X = xx-w:xx+w
                for Y = yy-w:yy+w
                    if Y > 0 && X > 0 && Isize(1) >= X && Isize(2) >= Y && Z > 0 && Isize(3) >= Z
                        binary_image(X,Y,Z) = 1;
                    end
                end
            end
        end
    end
    newim = double(binary_image);
