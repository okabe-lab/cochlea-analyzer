function [point3,C,point3_2] = ReObtainPointGroup2(point2,im0,xres,zres,figureon)
%*****************************************************************************
%This function uses "pointList", coordinates of points which may represent organ 
%of corti, "originalIm", image stack, and "PIXEL_WIDTH" and "Z_STEP", image resolution. The 
%outputs are "pointList", new coordinates of points, "coefficientsMat" result array of
%template matching, "pointList2", new coordinates of points with lax criteria.  
%*****************************************************************************

%% Make maximum projection image and detect reginal peaks
ppoint2 = point2*diag([1/xres,1/xres,1/zres]);
tempim = im0(:,:,round(min(ppoint2(:,3))):round(max(ppoint2(:,3))));
maxim = max(tempim,[],3);
maxim2 = medfilt2(maxim);
th3 = multithresh(maxim2(:),3);
thim = maxim2<th3(1);
temp2 = maxim2(thim>0);
meanBG = mean(temp2);
sdBG = std(temp2);
th4 = meanBG + sdBG*3;
thim2 = bwareaopen(maxim2>th4,2);
% ピーク検出
peakim2 = imregionalmax(maxim2);
peakim3 = peakim2.*thim2;

if figureon ==1
    figure
    imshow(peakim2,[])
end

im1max = maxim == 4095;

L0 = bwlabeln(im1max);
stats0 = regionprops(L0,'Area');
volume = cat(1,stats0.Area);
cell_width = 7/xres;

add_points = [];
for i = 1:size(volume,1)
    if volume(i) > 20*xres*xres*zres
        temp_im = L0 == i;
        temp_xyz = binaryToXYZ(temp_im); 
        temp1 = max(temp_xyz(:,1))-min(temp_xyz(:,1));
        temp2 = max(temp_xyz(:,2))-min(temp_xyz(:,2));
        if temp1 >= temp2 
            temp_num = ceil(temp1/cell_width);
            for j = 1:temp_num
                f1 = temp_xyz(:,1) >= min(temp_xyz(:,1))+(j-1)*cell_width;
                f2 = temp_xyz(:,1) < min(temp_xyz(:,1))+j*cell_width;
                f = f1.*f2;
                if size(temp_xyz(f>0,:),1)>1
                    add_points = [add_points; mean(temp_xyz(f>0,:))]; 
                else
                    add_points = [add_points; temp_xyz(f>0,:)];
                end
            end
        else
            temp_num = ceil(temp2/7);
            for j = 1:temp_num
                f1 = temp_xyz(:,2) >= min(temp_xyz(:,2))+(j-1)*cell_width;
                f2 = temp_xyz(:,2) < min(temp_xyz(:,2))+j*cell_width;
                f = f1.*f2;
                if size(temp_xyz(f>0,:),1)>1
                    add_points = [add_points; mean(temp_xyz(f>0,:))]; 
                else
                    add_points = [add_points; temp_xyz(f>0,:)];
                end
            end            
        end     
    end
end

if not(isempty(add_points))
    zlist2 = zeros(size(add_points,1),1);
    for i = 1:size(add_points,1)
        temp = round(add_points(i,1:2));
        temp2 = tempim(temp(1),temp(2),:);
        M = max(temp2);
        f = find(temp2 == M);
        zlist2(i,1) = f(end) + round(min(ppoint2(:,3))) -1;
    end
    add_points = [add_points zlist2];
end

%% Get coordinates of regional maxima
L = bwlabeln((peakim2)>0);
stats = regionprops(L,'centroid');
centroids = cat(1,stats.Centroid);

dth = 20; 
numth = 23; 
squareD = squareform(pdist(centroids));
dense = (sum(squareD<dth))';
f = (dense>numth);
g1 = centroids(f>0,:);
g2 = centroids(f==0,:);

dense1 = sum(repmat(f,1,size(f,1)).*(squareD<dth))';
dense2 = sum(repmat(f==0,1,size(f,1)).*(squareD<dth))';

g2_2 = centroids(dense1 < 5,:);  
g2_3 = centroids(dense2 > 12,:);
g3 = [g2_2; g2_3]; 

points = zeros(size(g3,1),2);
points(:,1) = g3(:,2);
points(:,2) = g3(:,1);

points = unique(round([points; ppoint2(:,1:2)]),'rows');

%% Get z coordinates
zlist = zeros(size(points,1),1);
for i = 1:size(points,1)
    temp = points(i,:);
    temp2 = tempim(temp(1),temp(2),:);
    [~,idx] = max(temp2);
    zlist(i,1) = idx + round(min(ppoint2(:,3))) -1;
end
point3 = [points zlist];

point3 = [point3; add_points];

[~,D] = knnsearch(ppoint2,point3);
point3(D > 100,:)=[];  

Y = pdist(point3,'euclid');
Z = linkage(Y,'single');
T = cluster(Z,'cutoff',24,'criterion','distance');
mem_num = zeros(max(T),1);
for i =1:max(T)
    mem_num(i,1) = sum(T==i);
end
f = mem_num > 20;
idx = find(f);
temp = [];
for i =1:sum(f)
    temp = [temp; point3(T==idx(i)>0,:)];
end
point3 = temp;

for ii = 1:5
coef = pca(point3);
point4 = point3*coef; 
zcenter = median(point4(:,3));
squareD = squareform(pdist(point4(:,1:2)));
tpoint3 = point3;
for i = 1:size(point4,1)
    temp = point4(i,:);
    temp2 = squareD(:,i);
    temp2(i,1) = inf;
    [minimum,idx] = min(temp2);
    if minimum < 3
        z1 = abs(temp(1,3)-zcenter);
        temp3 = point4(idx,:);
        z2 = abs(temp3(1,3)-zcenter);
        [~,idx2] = min([z1,z2]);
        if idx2 == 2
            tpoint3(i,:) = point3(idx,:);
        end
    end
end
point3 = unique(tpoint3,'rows');
end

    squareD = squareform(pdist(point3));
    tpoint3 = point3;
    for i = 1:size(point3,1)
        temp = round(point3(i,:));
        temp2 = squareD(:,i);
        temp2(temp2==0) = inf;
        [minimum,idx] = min(temp2);
        if minimum < 5
            intens1 = im0(temp(1),temp(2),temp(3));
            temp3 = round(point3(idx,:));
            intens2 = im0(temp3(1),temp3(2),temp3(3));
            [~,idx2] = max([intens1,intens2]);
            if idx2 == 2
                tpoint3(i,:) = temp3;
            end
        end
    end
    point3 = unique(tpoint3,'rows');


%% Template Matching
template = [166    166    175    175    175    143    150;...
    315    492    592    592    555    378    175;...
    525    687    945    945    776    718    412;...
    625    800    986   1005    986    800    669;...
    625    800    983   1012    966    800    669;...
    474    687    703    945    800    706    412;...
    250    380    528    582    582    328    234];
w = floor(size(template,1)/2);
C = normxcorr2(template,maxim2);
C = C(w+1:end-w,w+1:end-w);

noise_temp = [       60    61    59    52    58    55    53;...
    48    55    52    59    57    55    54;...
    64    62    59    69    53    56    58;...
    74    62    74   405   113    62    63;...
    69    56    77   103    73    65    75;...
    80    70    64    74    96    65    58;...
    68   100    85    82   136    95    65];

C2 = normxcorr2(noise_temp,maxim);
SE = strel('square',3);
C3 = imdilate(C2,SE);

corr = zeros(size(point3,2));
for i = 1:size(point3,1)
    corr(i,1) = C(round(point3(i,1)),round(point3(i,2)));
    corr(i,2) = C3(round(point3(i,1)),round(point3(i,2)));
end

if figureon == 1
    figure
    imshow(C,[])
    figure
    imshow(C3,[])
end

f = (corr(:,1)>0.1) + (corr(:,2)<0.6)>1;
point4 = point3(f>0,:);

Y = pdist(point4,'euclid');
Z = linkage(Y,'single');
T = cluster(Z,'cutoff',20,'criterion','distance');
mem_num = zeros(max(T),1);
for i =1:max(T)
    mem_num(i,1) = sum(T==i);
end
f = mem_num > 20;
idx = find(f);
temp = [];
for i =1:sum(f)
    temp = [temp; point4(T==idx(i)>0,:)];
end
point4 = temp;

if figureon == 1
    [size(point3) size(point4)]
    urata_newf46(point3);
    hold on
    scatter3(point4(:,1),point4(:,2),point4(:,3))
    hold off
    removed_noise1 = sum(corr(:,1)<=0.1)
    removed_noise2 = sum(corr(:,2)>=0.4)
end

point3_2 = point3; 
point3 = point4; 

point3 = point3*diag([xres, xres, zres]);
point3_2 = point3_2*diag([xres, xres, zres]);
