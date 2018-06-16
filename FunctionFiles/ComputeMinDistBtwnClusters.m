function [dist,X,Y,Z] = ComputeMinDistBtwnClusters(xyzlist1,xyzlist2,xyzlist3)
%*****************************************************************************
%This function uses "xyzlist1", "xyzlist2" and "xyzlist3", set of xyz 
%coordinates. The outputs are "dist", array of minimum distances between 
%clusters, "X", "Y" and "Z", array of componets of 3D vectors that 
%represents positional relation between clusters.  
%*****************************************************************************

%% Calculate minimum distance and positional relationship between clusters
% (Each cluster contains three points) 
xyzlist = zeros(size(xyzlist1));
num = size(xyzlist1,1);
for j = 1:num
    count = (j-1)*3+1;
    xyzlist(count,:) = xyzlist1(j,:);
    xyzlist(count+1,:) = xyzlist2(j,:);
    xyzlist(count+2,:) = xyzlist3(j,:);
end

D = squareform(pdist(xyzlist));
dist = zeros(num,num);
X = zeros(num,num);
Y = zeros(num,num);
Z = zeros(num,num);
for j = 1:num
    count1 = (j-1)*3+1;
    for k = 1:num
        count2 = (k-1)*3+1;
        tempD = D(count1:count1+2,count2:count2+2);
        [d,idx] = min(tempD(:));
        [row,col] = ind2sub([3 3],idx);
        vect = xyzlist(count2+col-1,:)-xyzlist(count1+row-1,:);
        dist(j,k) = d;
        X(j,k) = vect(1);
        Y(j,k) = vect(2);
        Z(j,k) = vect(3);
    end
end