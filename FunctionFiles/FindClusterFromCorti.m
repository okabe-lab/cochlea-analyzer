function idx = FindClusterFromCorti(candid_group,im0,th2,xres,zres)
%*****************************************************************************
%This function uses "candid_group", group of intensity peaks, "im0", image 
%stack, "th2", threshold for signal detection ,and "xres" and "zres", image
%resolution. The outputs are "idx", the most probable group representing
%organ of corti.
%*****************************************************************************

%% Measure inclination of approximation plane 
Isize = size(im0);
angles = zeros(size(candid_group,1),1);
for i = 1:size(candid_group,1)
    temp = pca(candid_group{i,1});
    angles(i,1) = dot(temp(:,3),[0 0 1]);
end

%% Measure expanse and density of peak points
latscat = zeros(size(candid_group,1),2);
f3 = zeros(size(candid_group,1),1);
f4 = zeros(size(candid_group,1),1);
for i = 1:size(candid_group,1)
    temp = pca(candid_group{i,1});
    temp2 = candid_group{i,1}*temp;
    f1 = temp2(:,1) > prctile(temp2(:,1),47);
    f2 = temp2(:,1) < prctile(temp2(:,1),53);
    f = f1.*f2;
    temp3 = temp2(f>0,:);
    latscat(i,1) = max(temp3(:,2))-min(temp3(:,2));
    latscat(i,2) = max(temp3(:,3))-min(temp3(:,3));
    f4(i,1) = max(temp2(:,1))-min(temp2(:,1));
    temp4 = (temp3 / temp)*diag([1/xres,1/xres,1/zres]);
    temp5 = round(median(temp4));    
    t1 = max(temp5(1)-40,1);
    t2 = min(temp5(1)+40,Isize(1));
    t3 = max(temp5(2)-40,1);
    t4 = min(temp5(2)+40,Isize(2));
    temp_im = zeros(Isize);
    temp_im(t1:t2,t3:t4,temp5(3)) = im0(t1:t2,t3:t4,temp5(3));
    temp_im = temp_im>th2*1.4;
    f3(i,:) = sum(temp_im(:)); 
end

%% Measure density of peak points
density_score = (f3>1000) + (f3>500);
f1 = (latscat(:,1) > 50).*(latscat(:,1) < 200);
f2 = (latscat(:,2) < 70);
latscat_score1 = f1.*f2;
latscat_score2 = (f4 > 350) - (f4 < 200);
angles_score = (angles > 0.9) + (angles > 0.8);

%% Template matching by image of hair cell
template = [166    166    175    175    175    143    150;...
    315    492    592    592    555    378    175;...
    525    687    945    945    776    718    412;...
    625    800    986   1005    986    800    669;...
    625    800    983   1012    966    800    669;...
    474    687    703    945    800    706    412;...
    250    380    528    582    582    328    234];
w = floor(size(template,1)/2);
corr = zeros(size(candid_group));
for i = 1:size(candid_group,1)
    ppoint2 = round(candid_group{i,1}*diag([1/xres 1/xres 1/zres]));
    tempim = im0(:,:,round(min(ppoint2(:,3))):round(max(ppoint2(:,3))));
    maxim = max(tempim,[],3);
    maxim2 = medfilt2(maxim);
    C = normxcorr2(template,maxim2);
    C = C(w+1:end-w,w+1:end-w);
    
    ind = sub2ind(size(C),ppoint2(:,1),ppoint2(:,2));
    tempcorr = C(ind);
    corr(i,1) = median(tempcorr);
end
corr_score = (corr>0.35) + (corr>0.4);

%% Score evaluation
score = corr_score + latscat_score1 + angles_score + density_score + latscat_score2;
[~,idx] = max(score);

