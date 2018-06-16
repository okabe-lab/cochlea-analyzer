function [occupied_im, centers_im, gap_coord, gap_centers, gap_region_props] = ...
    DetectCellLoss(out_cell_coord, ANALYZE_BEGIN, ANALYZE_END, IM_SIZE, EXCLUDE_RANGE)
%*****************************************************************************
%This function uses "out_cell_coord", cartesian coordinates of outer hair 
%cells in linearized image, "ANALYZE_BEGIN" and "ANALYZE_END", analyzing 
%range in y coordinate, "IM_SIZE" image size, and "EXCLUDE_RANGE", exclusion
%range of analysis. The outputs are "occupied_im", occupied area by outer 
%hair cells, "centers_im", centers of outer hair cells, "gap_coord" and 
%"gap_centers", location of estimated cell loss, "region_props", properties
%of void spaces.
%*****************************************************************************

%% Standerdize locations of outer hair cells
within = [ANALYZE_BEGIN ANALYZE_END];
xyz = out_cell_coord(:,1:2);
count = 0;
window = 4;

minmidmax = zeros(ANALYZE_END-ANALYZE_BEGIN+1,4);

for j = ANALYZE_BEGIN:ANALYZE_END
    if sum(ismember(EXCLUDE_RANGE,j))>0
        continue
    end
    
    f = (xyz(:,2)>j-window).*(xyz(:,2)<j+window);
    if size(xyz(f>0,1),1)<3
        continue
    end
    tempmin = min(xyz(f>0,1));
    tempmax = max((xyz(f>0,1)));
    if tempmax-tempmin < 10
        continue
    end
    count = count + 1;
    minmidmax(count,:) = [j,min(xyz(f>0,1)),(min(xyz(f>0,1))+max((xyz(f>0,1))))/2,max((xyz(f>0,1)))];
end
minmidmax(minmidmax(:,1)==0,:) = [];

w = 25;
num = floor((max(xyz(:,2))-min(xyz(:,2)))/w)-1;
stlist = min(xyz(:,2)):w:min(xyz(:,2))+w*num;
edlist = [min(xyz(:,2))+w:w:max(xyz(:,2))-w,max(xyz(:,2))+1];
mid1 = zeros(numel(stlist),2);
for j = 1:numel(stlist)
    f1 = minmidmax(:,1)>=stlist(j);
    f2 = minmidmax(:,1)< edlist(j);
    mid1(j,:) = [mean([stlist(j),edlist(j)]),mean(minmidmax(f1.*f2>0,3))];
end
mid1(isnan(mid1(:,2))>0,:) = [];
midline = interp1(mid1(:,1),mid1(:,2),(1:IM_SIZE(2))','linear','extrap');
midline2 = interp1(mid1(:,1),mid1(:,2),(1:IM_SIZE(2))');

temp2 = find(isnan(midline2)==0,1,'last');
midline(temp2+1:end,:) = midline(temp2);

xyz2 = xyz;
xyz2(:,1) = xyz2(:,1) - midline(round(xyz(:,2)),:);
temp = minmidmax(:,4)-minmidmax(:,2);

w = 25;
num = floor((max(xyz(:,2))-min(xyz(:,2)))/w)-1;
stlist = min(xyz(:,2)):w:min(xyz(:,2))+w*num;
edlist = [min(xyz(:,2))+w:w:max(xyz(:,2))-w,max(xyz(:,2))+1];
wid1 = zeros(numel(stlist),2);
for j = 1:numel(stlist)
    f1 = minmidmax(:,1)>=stlist(j);
    f2 = minmidmax(:,1)< edlist(j);
    wid1(j,:) = [mean([stlist(j),edlist(j)]),mean(temp(f1.*f2>0,1))];
end
wid1(isnan(wid1(:,2))>0,:) = [];
width = interp1(wid1(:,1),wid1(:,2),(1:IM_SIZE(2))');
temp1 = find(isnan(width)==0,1);
temp2 = find(isnan(width)==0,1,'last');
width(1:temp1-1,:) = width(temp1);
width(temp2+1:end,:) = width(temp2);

xyz2(:,1) = 10*xyz2(:,1)./width(round(xyz(:,2)),:);

int1 = zeros(size(xyz,1),2);
for j = 1:size(xyz,1)
    temp = xyz(j,:);
    f1 = abs(xyz(:,1)-temp(1))<4;
    f2 = abs(xyz(:,2)-temp(2))<10;
    f3 = abs(xyz(:,2)-temp(2))>4;
    temp2 = xyz((f1.*f2.*f3)>0,:);
    if ~isempty(temp2)
        int1(j,:) = [temp(2) mean(abs(temp2(:,2)-temp(2)))];
    end
end
int1(int1(:,1)==0,:)=[];

w = 100;
num = floor((max(xyz(:,2))-min(xyz(:,2)))/w)-1;
stlist = min(xyz(:,2)):w:min(xyz(:,2))+w*num;
edlist = [min(xyz(:,2))+w:w:max(xyz(:,2))-w,max(xyz(:,2))+1];
int2 = zeros(numel(stlist),2);
for j = 1:numel(stlist)
    f1 = int1(:,1)>=stlist(j);
    f2 = int1(:,1)< edlist(j);
    int2(j,:) = [mean([stlist(j),edlist(j)]),mean(int1(f1.*f2>0,2))];
end
int2(isnan(int2(:,2))>0,:) = [];
interval = interp1(int2(:,1),int2(:,2),(1:IM_SIZE(2))');
temp1 = find(isnan(interval)==0,1);
temp2 = find(isnan(interval)==0,1,'last');
interval(1:temp1-1,:) = interval(temp1);
interval(temp2+1:end,:) = interval(temp2);

templist = zeros(size(interval));
for j = 1:IM_SIZE(2)
    templist(j,1) = 5/interval(j);
end
translist = cumsum(templist);
xyz = xyz2;
xyz(:,2) = interp1((1:IM_SIZE(2))',translist,xyz(:,2));

%% Make standardized image indicating locations of outer hair cell
win = round(interp1(1:size(translist,1),translist,within,'linear','extrap'));
st2 = win(1);
ed2 = win(2);
simplane = zeros(17,round(ed2-st2));

if ~isempty(EXCLUDE_RANGE)
    for jj = 1:size(EXCLUDE_RANGE,1)
        tempom = EXCLUDE_RANGE(jj,:);
        temp = round(interp1(1:size(translist,1),translist,[min(tempom);max(tempom)])-st2);
        simplane(:,temp(1):temp(2)) = 1;
    end
end
simplane = padarray(simplane,[100,1000],1);
centers_im = simplane*255;
xyz2 = xyz+repmat([9 -st2],size(xyz,1),1);
for j = 1:size(xyz2,1)
    temp = round(xyz2(j,:)+[100 1000]);
    simplane(temp(1)-2:temp(1)+2,temp(2)-2:temp(2)+2)=1;
    centers_im(temp(1),temp(2)) = 255;
end
SE = strel('rectangle',[4,3]);
occupied_im = imclose(simplane,SE);
SE = strel('square',3);
occupied_im = imopen(occupied_im,SE);
occupied_im = occupied_im(101:end-100,1001:end-1000);
centers_im = centers_im(101:end-100,1001:end-1000);

%% Detect cell loss
temp = occupied_im(2:16,:)==0;
simplane4 = bwareaopen(temp,13);
emplist = GetCoordOfPositivePixels(simplane4);
emplist = emplist + repmat([0 st2],size(emplist,1),1);
gap_coord = emplist;
gap_coord(:,2) = interp1(translist,1:IM_SIZE(2),emplist(:,2));
gap_coord(:,1) = (gap_coord(:,1)-9).*width(round(gap_coord(:,2)),:)/10 + midline(round(gap_coord(:,2)) ...
    ,:) + 2;

CC = bwconncomp(simplane4,8);
gap_region_props = regionprops(CC,'Area','Centroid');
centerlist = cat(1,gap_region_props.Centroid);
centerlist = [centerlist(:,2) centerlist(:,1)];
gap_centers = centerlist + repmat([0 st2],size(centerlist,1),1);
gap_centers(:,2) = interp1(translist,1:IM_SIZE(2),gap_centers(:,2));
gap_centers(:,1) = (gap_centers(:,1)-9).*width(round(gap_centers(:,2)),:)/10 ...
    + midline(round(gap_centers(:,2)),:)+2;
