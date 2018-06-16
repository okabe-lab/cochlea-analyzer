function [comp1,vector] = MergeLines(centerline1_1,centerline2_1,vector1,vector2)
%*****************************************************************************
%This function uses "centerline1_1" and "centerline2_1", xyz coordinates of
%two lines to combine, and "vector1" and "vector2", vectors to combine 
%which indicate inclinations of apploximation planes corresponding to 
%"centerline1_1" and "centerline2_1", respectively. 
%*****************************************************************************

%Combine two fitting lines onto a single line and combine two set of 
%vectors indicating inclination of apploximation planes.

[idx,d] = knnsearch(centerline1_1,centerline2_1);
[~,idx2] = min(d);
idx3 = idx(idx2);
temp = centerline2_1(idx2,:) - centerline1_1(idx3,:);
temp2 = flip((0:1:size(centerline2_1,1)-idx2-1)/(size(centerline2_1,1)-idx2-1))';
temp3 = repmat(temp,size(centerline2_1,1)-idx2,1).*repmat(temp2,1,3);
comp1 = [centerline1_1(1:idx3,:); centerline2_1(idx2+1:end,:)-temp3];

temp = vector2(idx2,:) - vector1(idx3,:);

vector = [vector1(1:idx3,:); vector2(idx2+1:end,:)];
if idx3>2 && size(vector,1)> idx3+2
    vector(idx3-2,:) = vector(idx3-2,:) + temp/2/3;
    vector(idx3-1,:) = vector(idx3-1,:) + temp/2/2;
    vector(idx3,:) = vector(idx3,:) + temp/2;
    vector(idx3+1,:) = vector(idx3+1,:) - temp/2/2;
    vector(idx3+2,:) = vector(idx3+2,:) - temp/2/3;
    vector = vector./sqrt(sum(vector.^2,2));
end
