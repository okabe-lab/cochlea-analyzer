function [param,xyz4,norm2] = ObtainFeatureQuantities(xyz,lineindex,xyzlist2,idxlist ...
    ,L,scores2,leftend,outC)
%*****************************************************************************
%This function uses "xyz", coordinate of point, "lineindex", belonging row,
%"xyzlist2", list of xyz coordinates of peaks by template matching, 
%"idxlist" indice of points predicted to be outer hair cells, "L", label 
%matrix, "scores2", predition scores, "leftend", y coordinate of inner hair
%cell at the apical end, and "outC", matrix of correlation coefficients by 
%template matching. The outputs are "param", feature quantities of given 
%point, "xyz4", nearest peak of correlation coefficients from given point, 
%and "norm2", distance between "xyz" and "xyz4".
%*****************************************************************************

%% Extract feature quantity of given point
tempxyz = round(xyz);
idx3 = knnsearch(xyzlist2,xyz);
xyz4 = xyzlist2(idx3,:);
vector2 = xyz-xyz4;
norm2 = norm(vector2);
tempidx = L(xyz4(1),xyz4(2),xyz4(3));
tscore = scores2(idxlist == tempidx,2);

param = [xyz, xyz(2)-leftend, lineindex, outC(tempxyz(1),tempxyz(2),tempxyz(3))...
    vector2 norm2 tscore];