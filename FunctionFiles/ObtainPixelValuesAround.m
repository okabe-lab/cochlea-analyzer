function param = ObtainPixelValuesAround(xyz,compim2)
%*****************************************************************************
%This function uses "xyz", coordinates of point, "compim2", linearized
%image of organ of corti. The output is "param", reduced and normarized
%image for applying machine learning models.
%*****************************************************************************

%% Make reduced and normarized model centered on given coordinate
width_x = 12; width_y = 87;
width2 = 12; width3 = 8;

tcompim = padarray(compim2,[100 100 5]);
tempxyz = round(xyz);

stx = tempxyz(1)-width_x+100; edx = stx + width_x*2;
sty = tempxyz(2)-width_y+100; edy = sty + width_y*2;
stz = tempxyz(3)+5; edz = stz;
tempim = tcompim(stx:edx,sty:edy,stz:edz);
tempim(tempim==0) = mean(tempim(:));
tempim2 = histeq(double(tempim)/max(double(tempim(:))));

cent = [width_x+1,width_y+1];
tempim3 = imresize(tempim2(cent(1)-width2:cent(1)+width2,cent(2)-width2:cent(2)+width2) ...
    ,[width3 width3]);
param = tempim3(:)';

