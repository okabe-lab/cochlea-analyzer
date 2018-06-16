function pointList = ReObtainPointGroup(pointList,originalIm,PIXEL_WIDTH,Z_STEP)
%*****************************************************************************
%This function uses "pointList", coordinates of points which may represent organ 
%of corti, "originalIm", image stack, and "PIXEL_WIDTH" and "Z_STEP", image resolution. The 
%outputs are "pointList", new coordinates of points, "coefficientsMat" result array of
%template matching, "pointList2", new coordinates of points with lax criteria.  
%*****************************************************************************

%% Make maximum projection image and detect reginal peaks
ppoint2 = pointList*diag([1/PIXEL_WIDTH,1/PIXEL_WIDTH,1/Z_STEP]);
tempim = originalIm(:,:,round(min(ppoint2(:,3))):round(max(ppoint2(:,3))));
maxim = max(tempim,[],3);
maxim2 = medfilt2(maxim);
peakim2 = imregionalmax(maxim2);

%%  Make complementary points for saturated area in image stack
im1max = maxim == 4095;

L0 = bwlabeln(im1max);
stats0 = regionprops(L0,'Area');
volume = cat(1,stats0.Area);
cell_width = 7/PIXEL_WIDTH;

add_points = [];
for i = 1:size(volume,1)
    if volume(i) > 20*PIXEL_WIDTH*PIXEL_WIDTH*Z_STEP
        temp_im = L0 == i;
        temp_xyz = GetCoordOfPositivePixels(temp_im);
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

%% Remove noise
dth = 20;
numth = 23;
squareD = squareform(pdist(centroids));
dense = (sum(squareD<dth))';
f = (dense>numth);

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

%% Remove noise
point3 = RemoveSmallClusters(point3, 24, 21);
pointList = point3*diag([PIXEL_WIDTH, PIXEL_WIDTH, Z_STEP]);
