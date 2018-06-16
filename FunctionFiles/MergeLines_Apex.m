function [compline,newvect] = MergeLines_Apex(centerline1_1,centerline2_1,vector1,vector2)
%*****************************************************************************
%This function uses "centerline1_1" and "centerline2_1", xyz coordinates of
%two lines to combine, and "vector1" and "vector2", vectors to combine 
%which indicate inclinations of apploximation planes corresponding to 
%"centerline1_1" and "centerline2_1", respectively. The outputs are
%"compline", xyz coordinates of combined line, and "newvect", set of
%combined vectors. This function is used for image stack from apical end 
%which show rapid change of curvature.
%*****************************************************************************

%Combine two fitting lines onto a single line and combine two set of 
%vectors indicating inclination of apploximation planes.

temp1 = centerline1_1;
temp2 = centerline2_1;
temp3 = centerline1_1 + vector1;
temp4 = centerline2_1 + vector2;

[idx,d] = knnsearch(temp1,temp2);
[~,idx2] = min(d);
idx1 = idx(idx2);
temp1_2 = temp1(1:idx1,:);
temp3_2 = temp3(1:idx1,:);
temp2_2 = temp2(idx2+1:end,:);
temp4_2 = temp4(idx2+1:end,:);

temp1_3 = temp1_2(round(size(temp1_2,1)/2):end,:);
temp3_3 = temp3_2(round(size(temp1_2,1)/2):end,:);
temp2_3 = temp2_2(1:round(size(temp2_2,1)/2),:);
temp4_3 = temp4_2(1:round(size(temp2_2,1)/2),:);

coef = pca([temp1_3; temp2_3; temp3_3; temp4_3]);
ttemp1 = temp1_3*coef;
ttemp2 = temp2_3*coef;
ttemp3 = temp3_3*coef;
ttemp4 = temp4_3*coef;

temp12 = [ttemp1;ttemp2];
temp34 = [ttemp3;ttemp4];

qy1 = spline(temp12(1:80:end,1),temp12(1:80:end,2),temp12(:,1));
qz1 = spline(temp12(1:80:end,1),temp12(1:80:end,3),temp12(:,1));
qy2 = spline(temp34(1:80:end,1),temp34(1:80:end,2),temp34(:,1));
qz2 = spline(temp34(1:80:end,1),temp34(1:80:end,3),temp34(:,1));

comp1 = [temp12(:,1) qy1 qz1]/coef;
comp2 = [temp34(:,1) qy2 qz2]/coef;

comp1 = [temp1_2(1:round(size(temp1_2,1)/2)-1,:); comp1; temp2_2(round(size(temp2_2,1)/2)+1:end,:)];
comp2 = [temp3_2(1:round(size(temp1_2,1)/2)-1,:); comp2; temp4_2(round(size(temp2_2,1)/2)+1:end,:)];

lengs = sqrt(sum(diff(comp1).^2,2));
csleng = [0; cumsum(lengs)];

compline = interp1(csleng,comp1,0:1:csleng(end));
vectline = interp1(csleng,comp2,0:1:csleng(end));
newvect = movmean(vectline - compline,20);
newvect = newvect./sqrt(sum(newvect.^2,2));

vectline = compline + newvect;
lengs2 = sqrt(sum(diff(vectline).^2,2));
if min(lengs2) < 0.98
    [~,idx] = min(lengs2);
    st = max(1,idx-50);
    ed = min(idx+50,size(vectline,1));
    templine = vectline(st:ed,:);
    lengs3 = sqrt(sum(diff(templine).^2,2));
    csleng2 = [0; cumsum(lengs3)];
    templine2 = interp1(csleng2,templine,0:csleng2(end)/(size(templine,1)-1):csleng2(end));
    vectline(st:ed,:) = templine2;
    newvect = movmean(vectline - compline,10);
    newvect = newvect./sqrt(sum(newvect.^2,2));
end

%compline = compline/coef;
%newvect = newvect/coef;
