function inline = DetectInnerHairCells(compim,Mdl1,Mdl2)
%*****************************************************************************
%This function uses "compim", linearized image of organ of corti, and
%"Mdl1" and "Mdl2", machine learning models. The output is "inline", xyz
%coordinates of estimated locations of inner hair cells.
%*****************************************************************************

%% Template matching with typical image of inner hair cell
Isize = size(compim);
intemp1 = [428    541    854   1230   1307   1310   1253   1077    688    396    724
    582    661   1053   1457   1579   1787   1562   1495   1106    766   1126
    1048    991   1302   1583   1715   1789   1722   1620   1304   1458   1744
    1429   1328   1683   1990   2035   2127   2119   1750   1868   1816   1553
    2013   1639   2176   2700   2729   2511   2357   2182   2261   1846   1121
    2347   1542   1845   2292   2313   2381   2470   2410   1977   1597   1010
    1801    899   1328   1839   1918   1986   2000   2192   1724    997    821
    1011    544    773   1185   1050   1123   1077   1247    861    511    499
    721    374    385    450    500    630    462    412    280    263    300
    480    283    229    216    265    337    253    245    209    214    230
    247    205    197    195    212    229    214    200    228    206    220];

isize = size(intemp1);
inC = [];
inpeak = [];
for j = 1:Isize(3)
    temp = normxcorr2(intemp1,compim(:,:,j));
    inC = cat(3,inC,temp);
    
    peaktemp = imregionalmax(temp);
    inpeak = cat(3,inpeak,peaktemp);
end

inC = inC(round(isize(1)/2):round(isize(1)/2)+Isize(1)-1,round(isize(2)/2):round(isize(2)/2)+Isize(2)-1,:);
inpeak = inpeak(round(isize(1)/2):round(isize(1)/2)+Isize(1)-1,round(isize(2)/2)...
    :round(isize(2)/2)+Isize(2)-1,:);

if Isize(2) > 5000
    intemp2 =    [ 73    95   149   161   156   159   138   103    74
        95   133   170   167   162   157   134   120    99
        145   184   187   176   180   168   148   174   166
        201   255   272   267   237   247   225   209   234
        142   182   234   213   236   229   241   164   155
        90   108   106   118   108   114   115   107    93
        66    72    71    79    73    72    70    66    59
        53    58    62    65    61    56    52    59    57
        50    54    54    54    47    49    52    59    60];
    isize2 = size(intemp2);
    
    inC2 = [];
    inpeak2 = [];
    for j = 1:Isize(3)
        temp = normxcorr2(intemp2,compim(:,:,j));
        inC2 = cat(3,inC2,temp);
        
        peaktemp = imregionalmax(temp);
        inpeak2 = cat(3,inpeak2,peaktemp);
    end
    
    inpeak2 = inpeak2(round(isize2(1)/2):round(isize2(1)/2)+Isize(1)-1,round(isize2(2)/2) ...
        :round(isize2(2)/2)+Isize(2)-1,:);
    inpeak = [inpeak(:,1:5000,:) inpeak2(:,5001:end,:)];
end

inpeak3 = inpeak.*(inC>0);
inpeak3 = inpeak3.*(compim>0);
CC = bwconncomp(inpeak3>0,26);
L = labelmatrix(CC);

%% Feature quantity extraction
S = regionprops(CC,'Area','Centroid','PixelList','PixelIdxList');
param1 = [cat(1, S.Area), SwitchColumn1_2(cat(1, S.Centroid))];
gnum = double(max(L(:)));
param2 = zeros(gnum,9);
param3 = zeros(gnum,4);
for j = 1:gnum
    
    temp1 = SwitchColumn1_2(S(j).PixelList);
    temp2 = S(j).PixelIdxList;
    if size(temp1,1) > 1
        temp3 = inC(temp2).*double(compim(temp2));
        param2(j,1:3) = sum(temp1.*temp3./sum(temp3));
        
        [~,minidx] = min(temp1(:,3));
        [~,maxidx] = max(temp1(:,3));
    else
        param2(j,1:3) = temp1;
        minidx = 1; maxidx = 1;
    end
    param2(j,4:9) = [temp1(minidx,:) temp1(maxidx,:)];
    param3(j,1:2) = [max(inC(temp2)) double(max(compim(temp2)))];
end

leftend = 100;
rightend = Isize(2)-100;
interval = 150;
width = 150;
stlist = [(leftend:interval:rightend-width), rightend-width+1]';
c_ilist = zeros(size(stlist,1),4);
for j = 1:size(stlist,1)
    f = (param2(:,2) >= stlist(j)).*(param2(:,2) < stlist(j)+width);
    temp1 = param3(f>0,1:2);
    mxc = max(temp1(:,1));
    mxi = prctile(temp1(temp1(:,1)>mxc*0.8,2),90);
    mni = prctile(temp1(temp1(:,2)>0,2),5);
    c_ilist(j,:) = [stlist(j)+width/2 mxc mxi mni];
end

qc_ilist = interp1(c_ilist(:,1),c_ilist(:,2:4),(1:Isize(2))','linear','extrap');
ylist = round(param2(:,2));
param3(:,3) = param3(:,1)./qc_ilist(ylist,1);
param3(:,4) = (param3(:,2)-qc_ilist(ylist,3))./(qc_ilist(ylist,2)-qc_ilist(ylist,3));

interval = 50;
width = 50;
stlist = [(leftend:interval:rightend-width), rightend-width+1]';
bglist = zeros(size(stlist,1),2);
for j = 1:size(stlist,1)
    tempim = compim(round(Isize(1)/2):end,stlist(j):stlist(j)+width,:);
    tempim2 = tempim(:);
    tempim2(tempim2==0)=[];
    tempim3 = tempim2(tempim2<prctile(tempim2,30));
    bglist(j,:) = [stlist(j)+width/2 mean(tempim3)];
end

qbglist = interp1(bglist(:,1),bglist(:,2),(1:Isize(2))','linear','extrap');
param3_2 = param3(:,2)-qbglist(round(param2(:,2)));

temp = param2(:,1:3);
[IDX,D] = knnsearch(temp,temp,'K',2);
I = IDX(:,2);
param4 = [param1(I,:), param2(I,:) param3(I,:) D(:,2)];

intemp3 = [   2119   2667   3301   3459   2738   1994   2011   1219    731
    1763   1930   2465   3670   3167   2494   2338    861    538
    977    698   2584   3413   3507   1957   1513    363    356
    540    617   2992   3232   3338   2632   1567   1088    188
    271    629   2953   2863   3028   2699   2242   1565    256
    879    915   1798   2540   2934   2526   2767   1944    648
    1530   1878   2473   2607   2811   2553   2563   1821   1146
    1796   2296   2509   2667   2721   2678   2205   1652   1432
    1944   2088   2140   2392   2319   2469   2125   1527   1309
    1544   1769   1940   2032   2103   2346   1546    773   1006
    992   1505   1884   1857   2030   1782    508    406    611
    454    867    717   1134   1236    821    223    329    310
    461    452    356    263    281    381    173    171    200
    148    617    488    290    231    244    185    142    123
    106    346     88    127    160    129    150     88     85
    96     81     98    104    100    102     88     73     94
    98     90     88     77     96     94     92    104     92];

maxim = max(compim,[],3);
isize = size(intemp3);
Cinmx1 = normxcorr2(intemp3,maxim);
Cinmx2 = Cinmx1(round(isize(1)/2):end,round(isize(2)/2):end,:);
tempxy = round(param2(:,1:2))+repmat([2 0],size(param2,1),1);
param5 = Cinmx2(sub2ind(size(Cinmx2),tempxy(:,1),tempxy(:,2)));

params1 = [param1 param2 param3 param3_2 param4 param5];

%% First selection of candidate locations of inner hair cells
Mdl = Mdl1;
[~,score1] = predict(Mdl,params1);

width2 = 200;
num1 = round(Isize(2)/width2);
num2 = round(width2/3);
stlist = [0:floor(Isize(2)/num1):floor(Isize(2)/num1)*(num1-1) Isize(2)];
idxcell = cell(num1,1);
for j = 1:num1
    f = (params1(:,6) >= stlist(j)).*(params1(:,6) < stlist(j+1));
    tempscore = score1(:,2).*f;
    [~,tempidx] = sort(tempscore,'descend');
    idxcell{j,1} = tempidx(1:num2,:);
end
idx1_2 = cat(1,idxcell{:});
params2_1 = params1(idx1_2,:);
score1_2 = score1(idx1_2,2);

%% Additional feature quantity extraction
xyzlist1 = params2_1(:,5:7);
xyzlist2 = params2_1(:,8:10);
xyzlist3 = params2_1(:,11:13);
[dist,X,Y,Z] = ComputeMinDistBtwnClusters(xyzlist1,xyzlist2,xyzlist3);

param2_2 = zeros(size(params2_1,1),24);
paramidx = zeros(size(params2_1,1),7);

f0 = (dist > 4).*(dist <16);
for j = 1:size(params2_1,1)
    tempf = find(f0(:,j));
    tvect = [X(tempf,j), Y(tempf,j), Z(tempf,j)];
    normal = sqrt(sum(tvect.^2,2));
    nvect = tvect./normal;
    rad = acos(nvect(:,1));
    f = nvect(:,2)<0;
    rad(f>0) = 2*pi - rad(f>0);
    tempscore = score1_2(tempf);
    
    for k = 1:6
        f1 = rad >= (k-1)*(pi/3);
        f2 = rad < k*(pi/3);
        [M,idx] = max(tempscore.*f1.*f2);
        if M >0
            param2_2(j,(k-1)*4+1:k*4)=[tvect(idx,:),tempscore(idx)];
            paramidx(j,k) = tempf(idx);
        end
    end
end

[sortedD,I] = sort(dist,2);
IDX = I(:,2);
D = sortedD(:,2);
num = size(dist,1);
tempidx = sub2ind([num num],(1:num)',IDX);
vect = [X(tempidx) Y(tempidx), Z(tempidx)];
param2_3 = [vect score1_2(IDX) D];

params2 = [params2_1 param2_2 param2_3];

tcompim = padarray(compim,[100 100 5]);
iparams = zeros(size(params2,1),23*7);
for j = 1:size(params2,1)
    temp = round(params2(j,5:7));
    stx = temp(1)-40+100; edx = stx+69-1;
    sty = temp(2)-9+100; edy = sty+21-1;
    stz = temp(3)-5+5; edz = stz+5;
    tempim = histeq(max(tcompim(stx:edx,sty:edy,stz:edz),[],3));
    tempim = double(imresize(tempim,[23,7]));
    tempim = tempim/max(tempim(:));
    iparams(j,:) = tempim(:);
end

params3 = [params2 iparams];

%% Second selection of candidate locations of inner hair cells
Mdl = Mdl2;
[~,score] = predict(Mdl,params3);
score2 = score(:,2);

%% Detection of inner hair cells
% Apical end
Isize = size(compim);
th1 = 0.4;
th2 = 30;
num1 = 3;
th3 = 6;
th4 = 4;

points = params3(:,5:7);

f1 = points(:,2)<300;
f2 = score2 >= th1;
inxyz1 = points(f1.*f2>0,:);
score1 = score2(f1.*f2>0,:);

idx = knnsearch(inxyz1(:,1:2),[0 0]);
lu = inxyz1(idx,:);

dist = sort(squareform(pdist(inxyz1)));
dist = dist(2,:)';
inxyz1(dist>th2,:)=[];
score1(dist>th2,:)=[];

for jj = 1:3
    dist2 = squareform(pdist(inxyz1(:,1:2)));
    txyz = inxyz1;
    for j = 1:size(inxyz1,1)
        f = dist2(:,j);
        if sum(f<th4)>1
            [~,idx] = max(score1(f<th4,:));
            temp = inxyz1(f<th4,:);
            txyz(j,:) = temp(idx,:);
        end
    end
    inxyz1 = unique(txyz,'rows');
    inxyz1 = sortrows(inxyz1,2);
end

Y = pdist(inxyz1*diag([2,1,1]),'euclid'); 
Z = linkage(Y,'single'); 
T = cluster(Z,'cutoff',31,'criterion','distance');
mem_num = zeros(max(T),1);
for i =1:max(T)
    mem_num(i,1) = sum(T==i);
end
f = mem_num > 5;
idx = find(f);
candid_group = cell(sum(f),1);
xmin = zeros(sum(f),1);
for i =1:sum(f)
    temp = inxyz1(T==idx(i)>0,:);
    candid_group{i,1} = inxyz1(T==idx(i)>0,:);
    xmin(i) = min(temp(:,1));
end
[~,idx] = min(xmin);
inxyz1 = candid_group{idx,1};

idx = knnsearch(inxyz1(:,2),lu(1,2)+100);
temppoint = inxyz1(idx,:);
tempvect = temppoint(1,1:2)-lu(1,1:2);
slope = tempvect/tempvect(1,2);
tempxy = inxyz1(:,1:2)-repmat(lu(1,1:2),size(inxyz1,1),1);
difx = tempxy(:,1)-tempxy(:,2)*slope(1);
inxyz1(difx>15,:)=[];

if min(inxyz1(:,1)) > lu(1)
    inxyz1 = RemoveOutliersFromRow([inxyz1; lu],num1,th3);
    inxyz1 = unique([inxyz1; lu],'rows');
else
    inxyz1 = RemoveOutliersFromRow(inxyz1,num1,th3);
end
inxyz1 = sortrows(inxyz1,2);

f = sqrt(sum(diff(inxyz1(:,1:2)).^2,2))>25;
if sum(f)>0
    temp = find(f,1,'last');
    inxyz1 = inxyz1(temp+1:end,:);
end

if size(inxyz1,1)>=10
    inline1 = inxyz1(1:10,:);
else
    inline1 = inxyz1;
end

if min(inline1(:,2))>110
    load('machineLearningModels.mat','omdl1s');
    [params1,pred,~] = PredictOuterHairCells(compim(:,1:300,:),inline1,omdl1s{1,1});
    oxyz = params1(pred==1,5:7);
    f = oxyz(:,2)<min(inline1(:,2));
    if sum(f)>10
        temp = median(oxyz(f>0,:));
        temp = temp - [35 0 0];
        inline1 = [temp; inline1];
    end    
end

% Basal end
th1 = 0.4558;
th2 = 30;
num1 = 4;
th3 = 6;
th4 = 4;

f1 = points(:,2)>Isize(2)-200;
f2 = score2 >= th1;
inxyz2 = points(f1.*f2>0,:);
if size(inxyz2,1)<20
    temp = sort(score2(f1>0),'descend');
    th1 = temp(20,1);
    f2 = score2 >= th1;
    inxyz2 = points(f1.*f2>0,:);
end
inxyz0 = inxyz2;
score1 = score2(f1.*f2>0,:);

dist = sort(squareform(pdist(inxyz2)));
dist = dist(2,:)';
inxyz2(dist>th2,:)=[];

inxyz2(inxyz2(:,1)>median(inxyz0(:,1))+15,:) = [];
inxyz2(inxyz2(:,1)>Isize(1)/2+8,:) = [];

dist2 = squareform(pdist(inxyz2(:,1:2)));
txyz = inxyz2;
for j = 1:size(inxyz2,1)
    f = dist2(:,j);
    if sum(f<th4)>1
        [~,idx] = max(score2(f<th4,:));
        temp = inxyz2(f<th4,:);
        txyz(j,:) = temp(idx,:);
    end
end
inxyz2 = unique(txyz,'rows');

dist2_2 = squareform(pdist(inxyz2(:,2)));
th4_2 = 3;
txyz = inxyz2;
for j = 1:size(inxyz2,1)
    f = dist2_2(:,j);
    if sum(f<th4_2)>1
        [~,idx] = max(score2(f<th4_2,:));
        temp = inxyz2(f<th4_2,:);
        txyz(j,:) = temp(idx,:);
    end
end
inxyz2 = unique(txyz,'rows');
inxyz2 = sortrows(inxyz2,2);

inxyz2 = RemoveOutliersFromRow(inxyz2,num1,th3);
inline2 = inxyz2;

f = inxyz2(:,1) < prctile(inxyz2(:,1),80)-12;
inline2(f>0,:)=[];

% Middle part
interval = 200; 
th1 = 0.3715; 
th2 = 20; 
num1 = 4; 
th3 = 3.5;
th4 = 4.5;
xdist = 6;
ydist = 5;
gain = 1.5445;

leftend = inline1(1,:);
rightend = inline2(end,:);

stlist = [(leftend(1,2)+100:interval:rightend(1,2)-interval*2-101)...
    , rightend(1,2)-interval*2-100]';
inxyz4 = [];

for ii = 1:size(stlist,1)
    
    f1 = (points(:,2)>=stlist(ii)).*(points(:,2)<stlist(ii)+interval*2);
    f2 = score2 >= th1;
    inxyz3 = points(f1.*f2>0,:);
    score3 = score2(f1.*f2>0,:);

    if size(inxyz3,1) < round(interval/7)
        temp1 = points(f1>0,:);
        temp2 = score2(f1>0,:);
        [~,idx] = sort(temp2,'descend');
        tempnum = min(size(idx,1),round(gain*interval/7));
        inxyz3 = temp1(idx(1:tempnum),:);
        score3 = temp2(idx(1:tempnum),:);
    end

    dist = sort(squareform(pdist(inxyz3)));
    dist = dist(2,:)';
    inxyz3(dist>th2,:)=[];
    score3(dist>th2,:)=[];

    dist2 = squareform(pdist(inxyz3(:,1:2)));
    txyz = inxyz3;
    for j = 1:size(inxyz3,1)
        f = dist2(:,j);
        if sum(f<th4)>1
            [~,idx] = max(score3(f<th4,:));
            temp = inxyz3(f<th4,:);
            txyz(j,:) = temp(idx,:);
        end
    end
    inxyz3 = unique(txyz,'rows');
    inxyz3 = sortrows(inxyz3,2);
    inxyz4 = [inxyz4; inxyz3]; 
end

inxyz4 = unique(inxyz4,'rows');
IDX = knnsearch(points,inxyz4);
score4 = score(IDX,:);

dist3 = squareform(pdist(inxyz4(:,1:2)));
txyz = inxyz4;
for j = 1:size(inxyz4,1)
    f = dist3(:,j);
    if sum(f<th4)>1
        [~,idx] = max(score4(f<th4,:));
        temp = inxyz4(f<th4,:);
        txyz(j,:) = temp(idx(1),:);
    end
end
inxyz4 = unique(txyz,'rows');
inxyz4 = sortrows(inxyz4,2);

IDX = knnsearch(points,inxyz4);
score4 = score(IDX,:);
[~,idx] = max(score4);
stpoint = inxyz4(idx(1),:);
inxyz5 = DetectSingleRowOfPoints(stpoint,inxyz4,xdist,ydist);

inxyz5 = RemoveOutliersFromRow(inxyz5,num1,th3);
inline3 = inxyz5;

temp = [inline1; inline2; inline3];
inline = sortrows(unique(temp,'rows'),2);