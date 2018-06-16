function [params1,pred,score1,S,L,outC] = PredictOuterHairCells(compim2,inline2,Mdl)
%*****************************************************************************
%This function uses "compim2", linearized image of organ of corti,
%"inline2", xyz coordinates of inner hair cells, and "Mdl", machine
%learning model. The outputs are "params1", feature quantities of candidates
%for outer hair cells, "pred", prediction by the model, "score1", 
%prediction score calculated by the model, "S", properties of connected
%peak points, "L", label matrix of connected components, and "outC", matrix
%of correlation coefficients by template matching.
%*****************************************************************************

%% Template matching with typical image of outer hair cell
Isize = size(compim2);
outtemp1 = [    172    186    192    195    116     88     87
    330    543    533    487    380    242    135
    474    799    816    660    582    444    291
    610    999   1023    884    772    628    343
    565    960   1078    999    892    633    322
    435    676    850    892    754    495    297
    226    312    452    489    375    266    176];
isize = size(outtemp1);

outC = [];
outpeak = [];
for j = 1:Isize(3)
    temp = normxcorr2(outtemp1,compim2(:,:,j));
    outC = cat(3,outC,temp);
    
    peaktemp = imregionalmax(temp);
    outpeak = cat(3,outpeak,peaktemp);
end

outC = outC(round(isize(1)/2):round(isize(1)/2)+Isize(1)-1,round(isize(2)/2):round(isize(2)/2)...
    +Isize(2)-1,:);
outpeak = outpeak(round(isize(1)/2):round(isize(1)/2)+Isize(1)-1,round(isize(2)/2):round(isize(2)/2)...
    +Isize(2)-1,:);

leftend = round(min(inline2(:,2)));
rightend = round(max(inline2(:,2)));

outtemp2 = [     682   692   589   490   547   583   569   386   400   624
    315   425   621   761   725   747   573   425   248   411
    239   324   576   712   657   620   526   426   292   352
    370   265   520   611   643   635   614   595   428   310
    509   521   639   728   731   740   680   696   494   407];
isize2 = size(outtemp2);

outC2 = [];
outpeak3 = [];
for j = 1:Isize(3)
    temp = normxcorr2(outtemp2,compim2(:,:,j));
    outC2 = cat(3,outC2,temp);
    
    peaktemp = imregionalmax(temp);
    outpeak3 = cat(3,outpeak3,peaktemp);
end

outpeak3 = outpeak3(round(isize2(1)/2)-1:round(isize2(1)/2)+Isize(1)-2,round(isize2(2)/2) ...
    +1:round(isize2(2)/2)+Isize(2),:);
if (leftend + 4800) < Isize(2)
    outpeak = [outpeak(:,1:leftend+4800,:) outpeak3(:,leftend+4801:end,:)];
end

outpeak2 = outpeak.*(outC>0);
outpeak2 = outpeak2.*(compim2>0);
CC = bwconncomp(outpeak2>0,26);
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
        tempint = double(compim2(temp2));
        tempint = tempint-mean(tempint);
        tempint(tempint<0)=0;
        temp3 = outC(temp2).*tempint;
        if sum(temp3) <= 0
            temp3 = outC(temp2);
        end
        param2(j,1:3) = sum(temp1.*temp3./sum(temp3));        
        [~,minidx] = min(temp1(:,3));
        [~,maxidx] = max(temp1(:,3));
    else
        param2(j,1:3) = temp1;
        minidx = 1; maxidx = 1;
    end
    param2(j,4:9) = [temp1(minidx,:) temp1(maxidx,:)];
    param3(j,1:2) = [max(outC(temp2)) double(max(compim2(temp2)))];
end

interval = 150;
width = 150;
stlist = [(leftend:interval:rightend-width), rightend-width+1]';
if numel(stlist) == 1
    stlist = [stlist; stlist+10];
end
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
    tempim = compim2(round(Isize(1)/2):end,stlist(j):stlist(j)+width,:);
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
params1 = [param1 param2 param3 param3_2 param4];

%% First prediction on candidate locations of outer hair cells
[pred,score1] = predict(Mdl,params1);