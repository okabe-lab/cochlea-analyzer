function circle_center = ComputeCenterOfHelix(point_group,IS_FIRST)
%*****************************************************************************
%This function uses "point_group", cartesian coordinate of helically
%distributed points, "IS_FIRST", specify first round (1) or last round (0) 
%to calculate the center. The output is "circle_center", center coordinate 
%of fitted circle.
%*****************************************************************************

%% Extract one round of points from helically distributed points
if IS_FIRST == 1
    point_group = flipud(point_group);
end
squareD = squareform(pdist(point_group(:,1:2)));
distance1 = squareD(:,1);

movmean_dist = movmean(distance1,5);

interval = 500;
num = floor(max(point_group(:,4))/interval);
minvalue = zeros(num,1);
for i = 1:num
    f = (point_group(:,4)>(i-1)*interval).*(point_group(:,4)<i*interval);
    temp = movmean_dist(f>0,:);
    minvalue(i,1) = min(temp);
end

dif1_dist = [diff(minvalue)>0; 0];
dif2_dist = [0; diff(dif1_dist)==1];
min_point=find(dif2_dist,1);

if isempty(min_point)
    if IS_FIRST == 0
        first_circle = point_group(:,4) < min(point_group(:,4))+3000;
    else
        first_circle = point_group(:,4) > max(point_group(:,4))-3000;
    end    
else
    f = (point_group(:,4)>(min_point-1)*interval).*(point_group(:,4)<min_point*interval);
    temp = movmean_dist(f>0,:);
    idx = find(temp == minvalue(min_point,1));
    temp2 = point_group(f>0,:);
    if IS_FIRST == 0
        first_circle = point_group(:,4) < temp2(idx,4);
    else
        first_circle = point_group(:,4) > temp2(idx,4);
    end
end

X1 = point_group(first_circle>0,1:3);

%% Project point group onto principal component space
coeff = pca(X1);
X2 = X1*coeff;

%% Circle fitting
x = X2(:,1);
y = X2(:,2);
a0 = [mean(x),mean(y),max(x)-mean(x)];
f = @(a) norm((x-a(1)).^2 + (y-a(2)).^2- a(3).^2);
af1 = fminsearch(f, a0); 

center = [af1(1),af1(2),mean(X2(:,3))];
circle_center = center/coeff;