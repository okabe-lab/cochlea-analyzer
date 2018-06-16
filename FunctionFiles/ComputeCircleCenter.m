function center = ComputeCircleCenter(int_point, temp_center, number)
%*****************************************************************************
%This function uses "int_point", xy coordinates of points, "temp_center",
%initial value of center coordinate for circle fitting, "number", number of 
%partitions of points. The output is "center", center coordinates of fitted
%circle.
%*****************************************************************************

%% Extract representative points
temp = int_point(:,1:2)-repmat(temp_center,size(int_point,1),1);
rad0 = ComputeAngularCoord(temp);
dist0 = sqrt(temp(:,1).^2+temp(:,2).^2);

int_rad = min(rad0):(max(rad0)-min(rad0))/number:max(rad0);
rep_data = [];
for i = 1:size(int_rad,2)-1
    f1 = rad0>=int_rad(1,i);
    f2 = rad0< int_rad(1,i+1);
    f = f1.*f2;
    
    temp = prctile(dist0(f>0),90);
    f5_1 = dist0<temp+20;
    f5_2 = dist0>temp-50;
    f6 = f.*f5_1.*f5_2;
    if sum(f6) == 0
        continue
    end
    temp2 = prctile(dist0(f6>0),30);
    temp3 = knnsearch(dist0(f6>0),temp2);
    temp6 = int_point(f6>0,:);
    rep_data = [rep_data; temp6(temp3(1),:)];
end

%% Circle fitting
fun_1 = @(center) CircleFitRSS(rep_data,center);
center = fminsearch(fun_1,temp_center);


