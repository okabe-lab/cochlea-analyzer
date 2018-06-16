function compimage = DrawLinearizedIm(compimage,im0_1,o_xyz,o_xyz2,imf,Isize,difxyz)
%*****************************************************************************
%This function uses "compimage", linearized image under construction,
%"im0_1", original image stack, "o_xyz", xyz coordinates in linearized
%image, "o_xyz2", xyz coordinates in original image stack, "imf", array 
%indicating corresponding image stack for each pixel of linearized
%image, "Isize", size of original image stack, and "difxyz", shifts among
%image stacks. The output is "compimage", linearized image under
%construction.
%*****************************************************************************

im0_0 = padarray(im0_1,[1 1 1],'post');
for j = 1:size(o_xyz2,1)
    if imf(j) == 1
        temp = o_xyz2(j,:)-difxyz;
        if temp(1,1) >= 1 && temp(1,1) <= Isize(1) && temp(1,2) >= 1 && temp(1,2) ...
                <= Isize(2) && temp(1,3) >= 1 && temp(1,3) <= Isize(3)
            if temp(1) == round(temp(1)) && temp(2) == round(temp(2)) && temp(3) ...
                    == round(temp(3))
                compimage(o_xyz(j,1),o_xyz(j,2),o_xyz(j,3))=im0_0(round(temp(1,1)) ...
                    ,round(temp(1,2)),round(temp(1,3)));
            else
                temp2 = temp - floor(temp);
                temp3_1 = [(1-temp2(1))*(1-temp2(2))*(1-temp2(3)), (1-temp2(1))*temp2(2) ...
                    *(1-temp2(3));
                    temp2(1)*(1-temp2(2))*(1-temp2(3)), temp2(1)*temp2(2)*(1-temp2(3))];
                temp3_2 = [(1-temp2(1))*(1-temp2(2))*temp2(3), (1-temp2(1))*temp2(2)*temp2(3);
                    temp2(1)*(1-temp2(2))*temp2(3), temp2(1)*temp2(2)*temp2(3)];
                temp3 = cat(3,temp3_1,temp3_2);
                tempim = im0_0(floor(temp(1,1)):floor(temp(1,1))+1,floor(temp(1,2)):floor(temp(1,2)) ...
                    +1,floor(temp(1,3)):floor(temp(1,3))+1);
                temp4 = double(tempim).*temp3;
                compimage(o_xyz(j,1),o_xyz(j,2),o_xyz(j,3))=sum(temp4(:));
            end
        end
    end
end
