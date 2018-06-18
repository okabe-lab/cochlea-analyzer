%*****************************************************************************
%This script is continuation of "main1.m". This script perform detection of
%inner hair cells on 1st linearized image, create 2nd linearized image 
%based on inner hair cell row ("linearizedIm2.tif"). It also create an image with
%estimated inner hair cell locations ("innerHairCells.tif"), and excel file with
%their coordinates("innerHairCells.xlsx"). This script uses machine learning models
%("machineLearningModels.mat"), which we will provide on request.  
%*****************************************************************************

clearvars -except SubFolderPath
if not(exist('mainPath','var'))
    mainPath = fileparts(mfilename('fullpath')); % Get fullpath of this script
end
cd(mainPath)
load([SubFolderPath '\data.mat'],'linearizedIm1','shiftOfIms_pixel','PIXEL_WIDTH','Z_STEP' ...
    ,'centerLineOfCorti','inclineVectors','imFileNames','imSizeList','LINEAR_IM_WIDTH','LINEAR_IM_DEPTH')
load('machineLearningModels.mat','imdl1s','imdl2');
Isize = size(linearizedIm1);

%% Detection of inner hair cells from 1st image
inline = DetectInnerHairCells(linearizedIm1,imdl1s{1,1},imdl2);

inline3 = inline;
inline3(:,3) = movmean(inline3(:,3),4);
if inline3(end,2) < Isize(2)-45
    inline3 = [inline3; inline3(end,1) Isize(2) inline3(end,3)];
end
disp('Inner hair cells in 1st image detected...')

%% Create 2nd linearized image
vq = interp1(inline3(:,2),[inline3(:,1) inline3(:,3)],round(min(inline3(:,2)):max(inline3(:,2))) ...
    ,'linear','extrap');
qline = [round(min(inline3(:,2)):max(inline3(:,2)))' vq];
temp2 = qline - repmat([0,LINEAR_IM_WIDTH+1,LINEAR_IM_DEPTH+1],size(qline,1),1);
qline2 = centerLineOfCorti(temp2(:,1),:) - inclineVectors(temp2(:,1),:).*temp2(:,2);
qline2(:,3) = qline2(:,3) + temp2(:,3);

acdist = zeros(size(qline2,1),1);
for j = 2:size(qline2,1)
    acdist(j) = acdist(j-1) + norm(qline2(j,:) - qline2(j-1,:));
end

interval = 40;
xx = (0:interval:acdist(end)-interval)';
rep_point = zeros(size(xx,1),3);
rep_point2 = zeros(size(xx,1),3);
for j = 1:size(xx,1)
    idx = knnsearch(acdist,xx(j));
    rep_point(j,:) = qline2(idx,:); 
    rep_point2(j,:) = qline(idx,:); 
end
rep_point2 = SwitchColumn1_2(rep_point2);

rep_point2_3 = rep_point2;
rep_point2_3(:,1) = rep_point2(:,1)+30;

w = 40;
interval = 6;
rep_point2_2 = rep_point2;
for j = 1:size(rep_point2,1)
    temp = rep_point2(j,:);
    z = AdjustZCoord(linearizedIm1,temp,w,interval);    
    rep_point2_3(j,3) = z;    
end
rep_point2_3 = (rep_point2+rep_point2_3)/2;

centerline1 = interp1(rep_point2_3(:,2),[rep_point2_3(:,1) rep_point2_3(:,3)],1:size(linearizedIm1,2) ...
    ,'linear','extrap');
centerline1 = [(1:size(linearizedIm1,2))' centerline1];
temp2 = centerline1 - repmat([0,LINEAR_IM_WIDTH+1,LINEAR_IM_DEPTH+1],size(centerline1,1),1);
centerline2 = centerLineOfCorti(temp2(:,1),:) - inclineVectors(temp2(:,1),:).*temp2(:,2);
centerline2(:,3) = centerline2(:,3) + temp2(:,3);

acdist = zeros(size(centerline2,1),1);
for j = 2:size(centerline2,1)
    acdist(j) = acdist(j-1) + norm(centerline2(j,:) - centerline2(j-1,:));
end
interval = 1;
xx = (0:interval:acdist(end))';
centerline2_2 = zeros(size(xx,1),3);
centerline1_2 = zeros(size(xx,1),3);
for j = 1:size(xx,1)
    idx = find((acdist-xx(j))>0,1);
    tempdif = acdist(idx)-xx(j);
    if tempdif > 0
        centerline2_2(j,:) = centerline2(idx,:) + (centerline2(idx-1,:)-centerline2(idx,:))*tempdif;
    else 
        centerline2_2(j,:) = centerline2(end,:);
    end
    if tempdif > 0
        tempvect = (centerline1(idx-1,:)-centerline1(idx,:));
        tempvect = (tempvect/norm(tempvect))*tempdif;
        centerline1_2(j,:) = centerline1(idx,:) + tempvect;
    else 
        centerline1_2(j,:) = centerline2(end,:);
    end
end

inline_3 = interp1(rep_point2_2(:,2),[rep_point2_2(:,1) rep_point2_2(:,3)],centerline1_2(:,1) ...
    ,'linear','extrap');
inline_3 = [centerline1_2(:,1) inline_3];
temp = centerline1_2(1,:);
temp2 = inline_3(1:50,:)-temp;
temp3 = centerline1_2(1,:)-centerline1_2(50,:);
if temp3(2)<0
    [~,idx] = min(abs(dot(repmat(temp3,size(temp2,1),1),temp2,2)));
    tinline_3 = inline_3(idx:50,:);
    tinline_3 = interp1([1;50],[tinline_3(1,:); tinline_3(end,:)],1:50);
    inline_3(1:50,:) = tinline_3;
end

temp2 = inline_3 - repmat([0,LINEAR_IM_WIDTH+1,LINEAR_IM_DEPTH+1],size(inline_3,1),1);
temp3 = interp1(1:size(centerLineOfCorti,1),centerLineOfCorti,inline_3(:,1),'linear','extrap');
temp4 = interp1(1:size(inclineVectors,1),inclineVectors,inline_3(:,1),'linear','extrap');
inline_4 = temp3 - temp4.*temp2(:,2);
inline_4(:,3) = inline_4(:,3) + temp2(:,3);

num = 300;
acdist = zeros(num,1);
for j = 2:num
    acdist(j) = acdist(j-1) + norm(inline_4(j,:) - inline_4(j-1,:));
end
interval = acdist(end)/num;
xx = (0:interval:acdist(end))';
inline_5 = zeros(num,3);
for j = 1:num
    idx = find((acdist-xx(j))>0,1);
    tempdif = acdist(idx)-xx(j);
    if tempdif > 0
        inline_5(j,:) = inline_4(idx,:) + (inline_4(idx-1,:)-inline_4(idx,:))*tempdif;
    else 
        inline_5(j,:) = inline_4(num,:);
    end
end
inline_4(1:num,:) = inline_5;

num = 300;
acdist = zeros(num,1);
for j = 2:num
    acdist(j) = acdist(j-1) + norm(inline_4(num+j,:) - inline_4(num+j-1,:));
end
interval = acdist(end)/num;
xx = (0:interval:acdist(end))';
inline_5 = zeros(num,3);
for j = 1:num
    idx = find((acdist-xx(j))>0,1);
    tempdif = acdist(idx)-xx(j);
    if tempdif > 0
        inline_5(j,:) = inline_4(num+idx,:) + (inline_4(num+idx-1,:)-inline_4(num+idx,:)) ...
            *tempdif;
    else 
        inline_5(j,:) = inline_4(num,:);
    end
end
inline_4(num+1:num*2,:) = inline_5;
centerLineOfCorti2 = centerline2_2;
inline_2 = inline_4;

vector2 = inline_2-centerLineOfCorti2;
vector2 = vector2./sqrt(sum(vector2.^2,2));
inclineVectors2 = movmean(vector2,7);

zres2 = 1;
compimage = ones(size(centerLineOfCorti2,1),LINEAR_IM_WIDTH*2+1,LINEAR_IM_DEPTH*2+1);
o_xyz = GetCoordOfPositivePixels(compimage);
temp = o_xyz - repmat([0,LINEAR_IM_WIDTH+1,LINEAR_IM_DEPTH+1],size(o_xyz,1),1);
o_xyz2 = centerLineOfCorti2(temp(:,1),:) + inclineVectors2(temp(:,1),:).*temp(:,2);
o_xyz2(:,3) = o_xyz2(:,3) + temp(:,3)*zres2;
o_xyz3 = o_xyz2*diag([1/PIXEL_WIDTH,1/PIXEL_WIDTH,1/Z_STEP]);

imf = zeros(size(o_xyz2,1),size(imFileNames,1));
acf = zeros(size(o_xyz2,1),1);
for j = 1:size(imFileNames,1)
    f1 = (o_xyz3(:,1)-shiftOfIms_pixel(j,1) >= 1).*(o_xyz3(:,1)-shiftOfIms_pixel(j,1) <= imSizeList(j,1));
    f2 = (o_xyz3(:,2)-shiftOfIms_pixel(j,2) >= 1).*(o_xyz3(:,2)-shiftOfIms_pixel(j,2) <= imSizeList(j,2));
    f3 = (o_xyz3(:,3)-shiftOfIms_pixel(j,3) >= 1).*(o_xyz3(:,3)-shiftOfIms_pixel(j,3) <= imSizeList(j,3));
    f00 = f1.*f2.*f3;
    imf(:,j) = (f00 - acf)>0;
    acf = (acf + f00)>0;
end

compimage = zeros(size(centerLineOfCorti2,1),LINEAR_IM_WIDTH*2+1,LINEAR_IM_DEPTH*2+1);
for j = 1:size(imFileNames,1)
    load([SubFolderPath '\' TrimTif(imFileNames{j,1}) '.mat'],'-mat','processedIm');
    compimage = DrawLinearizedIm(compimage,processedIm,o_xyz,o_xyz3,imf(:,j),imSizeList(j,:) ...
        ,shiftOfIms_pixel(j,:));
end
linearizedIm2 = flipud(permute(uint16(compimage),[2 1 3]));

ImWrite3D(linearizedIm2,[SubFolderPath '\linearizedIm2.tif']);
disp('Re_linearization completed...')

%% Detection of inner hair cells in 2nd image
inline = DetectInnerHairCells(linearizedIm2,imdl1s{1,1},imdl2);

innerCells2 = inline;
innerCells2(:,3) = movmean(innerCells2(:,3),4);

incim = DrawMarkedIm(linearizedIm2,round(innerCells2),2,1);
incim = uint16(incim*1000);
disp('Inner hair cells in 2nd image detected!')

%% Save and export data
save([SubFolderPath '\data.mat'],'linearizedIm2','centerLineOfCorti2','inclineVectors2','innerCells2','-append')
xlswrite([SubFolderPath '\innerHairCells.xlsx'],sortrows(SwitchColumn1_2(round(innerCells2))),1)
ImWrite3D(incim,[SubFolderPath '\innerHairCells.tif']);
