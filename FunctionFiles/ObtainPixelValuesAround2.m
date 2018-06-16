function param = ObtainPixelValuesAround2(xyz,compim2)
%*****************************************************************************
%This function uses "xyz", coordinates of point, "compim2", linearized
%image of organ of corti. The output is "param", reduced and normarized
%image for applying machine learning models.
%*****************************************************************************

%% Make reduced and normarized model centered on given coordinate
tcompim = padarray(compim2,[100 100 5]);
wx = 69; wy = 39; wz = 1;

temp = round(xyz);
stx = temp(1)-floor(wx/2)+100; edx = stx+wx-1;
sty = temp(2)-floor(wy/2)+100; edy = sty+wy-1;
stz = temp(3)-floor(wz/2)+5; edz = stz+wz-1;
tempim = double(max(tcompim(stx:edx,sty:edy,stz:edz),[],3));
param = histeq(tempim/max(tempim(:)));

