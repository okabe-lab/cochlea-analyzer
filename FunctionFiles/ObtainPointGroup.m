function [points,th2,add_points] = ObtainPointGroup(im,xres,zres)
%*****************************************************************************
%This function uses "im", image stack, "xres" and "zres", image resolution.
%The outputs are "points", coordinates of intensity peaks, "th2", threshold
%for signal detection, and "add_points", complementary points for saturated 
%area in image stack.
%*****************************************************************************

%% Preprocess image stack and detect reginal maxima
im = medfilt3(im);
th1 = multithresh(im(:),2);
temp = im<th1(1);
temp2 = im(temp>0);
meanBG = mean(temp2);
sdBG = std(temp2);
th2 = meanBG + sdBG*5;
im2 = bwareaopen(im>th2,30);
peakim0 = imregionalmax(im);
peakim = peakim0.*im2;

%% Make complementary points for saturated area in image stack
im1max = im == 4095;

L0 = bwlabeln(im1max);
stats0 = regionprops(L0,'Area');
volume = cat(1,stats0.Area);
cell_width = 7/xres;
add_points = [];
for i = 1:size(volume,1)
    if volume(i) > 1000*xres*xres*zres
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

%% Get coordinates of regional maxima
L = bwlabeln((peakim)>0);
stats = regionprops(L,'centroid');
centroids = cat(1,stats.Centroid);

if size(centroids,1) > 20000 
 im2 = bwareaopen(im>th2*1.4,100);
 peakim = peakim0.*im2;
 L = bwlabeln((peakim)>0);
 stats = regionprops(L,'centroid');
 centroid2 = cat(1,stats.Centroid);
else
 centroid2 = centroids;
end

if size(centroid2,1) > 20000 
 im2 = bwareaopen(im>th2*2,100);
 peakim = peakim0.*im2;
 L = bwlabeln((peakim)>0);
 stats = regionprops(L,'centroid');
 centroid2 = cat(1,stats.Centroid);
end

%% Remove noise
temp = (sum(sort(squareform(pdist(centroid2)))<20))';
f = (temp<3) + (temp>29);
centroid2(f>0,:) =[];

points = zeros(size(centroid2,1),3);
points(:,1) = centroid2(:,2);
points(:,2) = centroid2(:,1);
points(:,3) = centroid2(:,3);

%% Add complemantary points 
points = [points; add_points];
% Convert to actual distance
points = points * diag([xres,xres,zres]);