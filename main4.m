%*****************************************************************************
%This script is continuation of "main3.m". The script analyzes degree of 
%cell loss. The script get coordinates of inner and outer hair cells from 
%excel files "innerHairCells.xlsx" and "outerHairCells.xlsx", respectively. These excel files
%are made by the preceding script "main3.m", and they can be modified by
%manual inspection before running this script. This script uses template
%for 3D registration ("cochlearTemplate.mat"), which we will provide on request. 
%The result of anlysis will save in MAT-file "result2.mat". The script also
%makes "standardized_OuterHairCells.tif", a normalized image indicating locations of outer
%hair cells.
%*****************************************************************************

%% Information input
% Specify analysis range by y coordinates of linearized image
analysis_range = []; % for example [33 4860]
omit = []; % for example [3247 3330]

%% Get coordinates of inner hair cells
load([SubFolderPath '\data.mat'],'linearizedIm2','LINEAR_IM_WIDTH','LINEAR_IM_DEPTH','centerLineOfCorti2','inclineVectors2')
Isize = size(linearizedIm2);
tempxyz = xlsread(fullfile(SubFolderPath,'innerHairCells.xlsx'));
tempxyz = tempxyz(:,1:3);
inline1 = sortrows(SwitchColumn1_2(tempxyz),2);

%% Transform coordinates in linearized image into those of original image stack
inst = inline1(1,:);
ined = inline1(end,:);
interval = 50;

stlist = round(inst(1,2)+interval/2:interval:ined(1,2))';
inline2 = zeros(numel(stlist),3);
for i = 1:numel(stlist)
    tempy = stlist(i);
    f = (inline1(:,2) > tempy-interval/2).*(inline1(:,2) < tempy+interval/2);
    if sum(f)>0
        tempxyz = inline1(f>0,:);
        inline2(i,:) = mean(tempxyz);
    end
end
f = inline2(:,1)==0;
inline2(f>0,:) = [];
inline2 = [inst; inline2; ined];

inline3 = interp1(inline2(:,2),inline2,(inst(1,2):ined(1,2))');
qline = [inline3(:,2) inline3(:,1) inline3(:,3)];
temp2 = qline - repmat([0,LINEAR_IM_WIDTH+1,LINEAR_IM_DEPTH+1],size(qline,1),1);
qline2 = centerLineOfCorti2(round(temp2(:,1)),:) - inclineVectors2(round(temp2(:,1)),:).*temp2(:,2);
qline2(:,3) = qline2(:,3) + temp2(:,3);

acdist = zeros(size(qline2,1),1);
for j = 2:size(qline2,1)
    acdist(j) = acdist(j-1) + norm(qline2(j,:) - qline2(j-1,:));
end

inlength = acdist(end);

indist = interp1(inst(1,2):ined(1,2),acdist,1:Isize(2),'linear','extrap')';
idist = interp1(inst(1,2):ined(1,2),acdist,inline1(:,2),'linear','extrap');

%% Set cylindrical coordinate system along modiolus
interval = 50;
xx2 = (0:interval:acdist(end))';
rep_point = zeros(size(xx2,1),3);
rep_point2 = zeros(size(xx2,1),3);
for i = 1:size(xx2,1)
    idx = knnsearch(acdist,xx2(i));
    rep_point(i,:) = qline2(idx,:);
end

center = ComputeCenterOfHelix([rep_point xx2],0);
center2 = ComputeCenterOfHelix([rep_point xx2],1);
centerdir = center-center2;
endpts = [center+centerdir*1.5;center2-centerdir*1.5];

a = cross(centerdir/norm(centerdir),[0,0,1]);
deg = acos(dot(centerdir/norm(centerdir),[0,0,1]));
xyz2 = RodriguesRotation(rep_point(:,1:3),a,deg);

ff = @(a) ComputeRssOfCochlearHelixFit(xyz2,a(1),a(2));
aa = fminsearch(ff,[0,0]);

temp = RodriguesRotation(xyz2,[1,0,0],aa(1));
xyz3 = RodriguesRotation(temp,[0,1,0],aa(2));

[~,cent1,cent2,sita,dist] =ComputeRssOfCochlearHelixFit(xyz2,aa(1),aa(2));
xyz4 = [1 0 0 -cent1; 0 1 0 -cent2; 0 0 1 -xyz3(1,3); 0 0 0 1]*[xyz3 ones(size(xyz3,1),1)]';
xyz4 = xyz4(1:3,:)';
if xyz4(end,3)<0
    xyz4(:,3)=-xyz4(:,3);
    xyz4(:,1)=-xyz4(:,1);
end
qline4 = [dist.*cos(sita) dist.*sin(sita) xyz4(:,3)];

temp2 = round(SwitchColumn1_2(inline1)) - repmat([0,LINEAR_IM_WIDTH+1,LINEAR_IM_DEPTH+1],size(inline1,1),1);
qline = centerLineOfCorti2(temp2(:,1),:) - inclineVectors2(temp2(:,1),:).*temp2(:,2);
qline(:,3) = qline(:,3) + temp2(:,3);
temp1 = RodriguesRotation(qline,a,deg);
temp2 = RodriguesRotation(temp1,[1,0,0],aa(1));
temp3 = RodriguesRotation(temp2,[0,1,0],aa(2));
temp4 = [1 0 0 -cent1; 0 1 0 -cent2; 0 0 1 -xyz3(1,3); 0 0 0 1]*[temp3 ones(size(temp3,1),1)]';
inline3 = temp4(1:3,:)';
if inline3(end,3)<0
    inline3(:,3)=-inline3(:,3);
    inline3(:,1)=-inline3(:,1);
end

sita1 = zeros(size(inline3,1),1);
dist1 = zeros(size(inline3,1),1);
for i = 1:size(inline3,1)
    temp = inline3(i,1:2);
    dist1(i,1) = norm(temp);
    
    temp = temp/norm(temp);
    
    if i == 1
        v0 = temp;
        v00 = temp;
        sita1(i,1) = 0;
        judge = cross([v00 0],[temp 0]);
        continue
    end
    
    temp2 = acos(dot(v0,temp));
    
    v0 = temp;
    sita1(i,1) = sita1(i-1,1) + temp2;
end
inline4 = [dist1.*cos(sita1) dist1.*sin(sita1) inline3(:,3)];
isita = interp1(inline1(:,2),sita1,(inst(1,2):ined(1,2))');

%% Get coordinates of inner hair cells
xyz = xlsread([SubFolderPath '\outerHairCells.xlsx']);
outline1 = sortrows(SwitchColumn1_2(xyz),2);

tlist = outline1;
temp2 = tlist - repmat([LINEAR_IM_WIDTH+1,0,LINEAR_IM_DEPTH+1],size(tlist,1),1);
tlist2 = centerLineOfCorti2(round(temp2(:,2)),:) - inclineVectors2(round(temp2(:,2)),:).*temp2(:,1);
tlist2(:,3) = tlist2(:,3) + temp2(:,3);

temp1 = RodriguesRotation(tlist2,a,deg);
temp2 = RodriguesRotation(temp1,[1,0,0],aa(1));
temp3 = RodriguesRotation(temp2,[0,1,0],aa(2));
temp4 = [1 0 0 -cent1; 0 1 0 -cent2; 0 0 1 -xyz3(1,3); 0 0 0 1]*[temp3 ones(size(temp3,1),1)]';
outline2 = temp4(1:3,:)';
if outline2(end,3)<0
    outline2(:,3)=-outline2(:,3);
    outline2(:,1)=-outline2(:,1);
end

odist = interp1(inst(1,2):ined(1,2),acdist,outline1(:,2),'linear','extrap');

%% 3D Registration with template data
load('cochlearTemplate.mat','z_temp','s_temp')
sita0 = (0:0.01:pi*4)';
midsita = max(s_temp)/2;

temp = interp1(sita1,inline4(:,3),sita0);
temp(isnan(temp),:) = [];
tz = temp;

t1 = z_temp(100:min(numel(tz),numel(z_temp))-50,:);
t2 = tz;
r = normxcorr2(t1,t2);
[~,I] = max(r);
lag = I - numel(t1) - 99;
csita1 = sita1 - lag*0.01;

temp1 = NaN(numel(s_temp),1);
st1 = max(1,lag);
ed1 = min(numel(s_temp)+lag-1,numel(tz));
st2 = max(1,-lag);
temp1(st2:st2+(ed1-st1)) = tz(st1:ed1);

temp = z_temp-temp1;
f = @(x) nansum((temp - x).^2);
zlag = fminsearch(f,0);
cz = inline4(:,3) + zlag;

midy = interp1(csita1,idist,midsita);

%% Analyze cell loss
if ~isempty(omit)
    f = zeros(size(odist));
    for k = 1:size(omit,1)
        com = interp1(inline1(:,2),idist,omit(k,:));
        f((odist>com(1)).*(odist<=com(2))>0)=1;
    end
    odist(f>0) = [];
end

if isempty(analysis_range)
    st = round(min(inline1(:,2)));
    ed = round(max(outline1(:,2)));
else
    st = min(analysis_range);
    ed = max(analysis_range);
end

[occupied_im, simplane2, emplist2, centerlist, S] = DetectCellLoss(outline1, st, ed, Isize, omit);
emptydist = interp1(inline1(:,2),idist,centerlist(:,2));
emptylist = cat(1, S.Area);
enum = round(emptylist/25);
om0 = omit/max(odist);

%% Measure distance with basal end as starting point
range_c = 3200;
idist2 = -(idist - midy)+range_c;
odist2 = -(odist - midy)+range_c;
edist2 = interp1(inline1(:,2),idist2,emplist2(:,2));
rend = interp1(inline1(:,2),idist2,ed,'linear','extrap');
edist3 = interp1(inline1(:,2),idist2,centerlist(:,2));

if ~isempty(omit)
    com = omit;
    fi = zeros(size(idist2));
    fe = zeros(size(edist2));
    fo = zeros(size(odist2));
    for k = 1:size(omit,1)
        com(k,:) = interp1(inline1(:,2),idist2,omit(k,[2 1]));
        fi((idist2>com(k,1)).*(idist2<=com(k,2))>0)=1;
        fe((edist2>com(k,1)).*(edist2<=com(k,2))>0)=1;
        fo((odist2>com(k,1)).*(odist2<=com(k,2))>0)=1;
    end
    idist2(fi>0) = [];
    edist2(fe>0) = [];
    odist2(fo>0) = [];
else
    com = [];
end

idist2(idist2<rend)=[];
odist2(odist2<rend)=[];
edist2(edist2<rend)=[];
enum2 = enum(edist3>=rend);

%% Evaluate degrees of cell loss in equally divided segments 
interval = 50;
xx = 0:interval:6000;

num = NaN(numel(xx)-1,1);
for i = 1:numel(num)
    tst = xx(i);
    ted = xx(i+1);
    f1 = idist2 >= tst;
    f2 = idist2 < ted;
    tin = idist2(f1.*f2>0);
    f3 = odist2 >= tst;
    f4 = odist2 < ted;
    tout = odist2(f3.*f4>0);
    tt = [tin; tout];
    if isempty(tt)
        continue
    elseif (max(tt)-min(tt)) < interval*0.5
        continue
    end
    f5 = edist2 >= tst;
    f6 = edist2 < ted;
    tem = edist2(f5.*f6>0);
    num(i) = round(numel(tem)/25);
end
void50 = num;

%% Calculate overall loss ratio 
conum = numel(odist2);
Ratio2 = sum(enum2)/(conum+sum(enum2));

%% Save and export results
tempim = cat(3,occupied_im,simplane2,zeros(size(occupied_im)));
imwrite(tempim(2:16,:,:),[SubFolderPath '\standardized_OuterHairCells.tif'])
save([SubFolderPath '\analyzeResults.mat'],'void50','Ratio2','enum2');
