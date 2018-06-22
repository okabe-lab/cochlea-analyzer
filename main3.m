%*****************************************************************************
%This script is continuation of "main2.m". This script perform detection of
%outer hair cells on 2nd linearized image. It creat an image with estimated
%outer hair cell locations ("outerHairCells.tif"), excel files with estimated
%locations of inner and outer hair cells ("outerHairCells.xlsx"). This script uses 
%machine learning models ("machineLearningModels.mat").
%*****************************************************************************

clearvars -except SubFolderPath
if not(exist('mainPath','var'))
    mainPath = fileparts(mfilename('fullpath')); % Get fullpath of this script
end
cd(mainPath)
load([SubFolderPath '\data.mat'],'linearizedIm2','innerCells2')
load('machineLearningModels.mat','omdl1s','omdl2','net1','net2');

%% First selection of candidate locations of outer hair cells
Mdl = omdl1s{1,1};
[params1,pred,score1,S, L, outC] = PredictOuterHairCells(linearizedIm2,innerCells2,Mdl);
disp('Cell candidates obtained...')

Isize = size(linearizedIm2);
width = 200;
num1 = round(Isize(2)/width);
num2 = round(width/2);
stlist = [0:floor(Isize(2)/num1):floor(Isize(2)/num1)*(num1-1) Isize(2)];
idxcell = cell(num1,1);
for j = 1:num1
    f = (params1(:,6) >= stlist(j)).*(params1(:,6) < stlist(j+1));
    tempscore = score1(:,2).*f;
    [~,tempidx] = sort(tempscore,'descend');
    idxcell{j,1} = tempidx(1:num2,:);
end
idxlist = cat(1,idxcell{:});

%% Additional feature quantity extraction
params2 = params1(idxlist,:);
score2 = score1(idxlist,2);

xyzlist1 = params2(:,5:7);
xyzlist2 = params2(:,8:10);
xyzlist3 = params2(:,11:13);
[dist,X,Y,Z] = ComputeMinDistBtwnClusters(xyzlist1,xyzlist2,xyzlist3);

param5 = zeros(size(params2,1),24);
paramidx = zeros(size(params2,1),7);

f0 = (dist > 4).*(dist <16);
for j = 1:size(params2,1)
    tempf = find(f0(:,j));
    tvect = [X(tempf,j), Y(tempf,j), Z(tempf,j)];
    normal = sqrt(sum(tvect.^2,2));
    nvect = tvect./normal;
    rad = acos(nvect(:,1));
    f = nvect(:,2)<0;
    rad(f>0) = 2*pi - rad(f>0);
    tempscore = score2(tempf);
    
    for k = 1:6
        f1 = rad >= (k-1)*(pi/3);
        f2 = rad < k*(pi/3);
        [M,idx] = max(tempscore.*f1.*f2);
        if M >0
            param5(j,(k-1)*4+1:k*4)=[tvect(idx,:),tempscore(idx)];
            paramidx(j,k) = tempf(idx);
        end
    end
end
[sortedD,I] = sort(dist,2);
IDX = I(:,2);
D = sortedD(:,2);
num = size(dist,1);
tempidx = sub2ind([num num],(1:num)',IDX);
vect = [X(tempidx) Y(tempidx), Z(tempidx)];
param6 = [vect score2(IDX) D];
paramidx(:,7) = IDX;

tcompim = padarray(linearizedIm2,[100 100 5]);
param7 = zeros(size(params2,1),21*7);
for j = 1:size(params2,1)
    temp = round(params2(j,5:7));
    stx = temp(1)-28+100; edx = stx+63-1;
    sty = temp(2)-10+100; edy = sty+21-1;
    stz = temp(3)-1+5; edz = stz+3-1;
    tempim = histeq(max(tcompim(stx:edx,sty:edy,stz:edz),[],3));
    tempim = imresize(tempim,[21,7]);
    tempim = double(tempim);
    param7(j,:) = tempim(:)/prctile(tempim(:),90);
end

params2_2 = [params2 param5 param6 param7];

xyzlist = params2_2(:,5:7);
for j = 1:size(idxlist)
    temp = SwitchColumn1_2(S(idxlist(j)).PixelList);
    [idx,~] = knnsearch(temp,xyzlist(j,:));
    xyzlist(j,:) = temp(idx,:);
end

tcompim = padarray(linearizedIm2,[100 100 5]);
wx = 69; wy = 39; wz = 3;
iparams = zeros(69,39,1,size(xyzlist,1));
for j = 1:size(xyzlist,1)
    temp = round(xyzlist(j,:));
    stx = temp(1)-floor(wx/2)+100; edx = stx+wx-1;
    sty = temp(2)-floor(wy/2)+100; edy = sty+wy-1;
    stz = temp(3)-floor(wz/2)+5; edz = stz+wz-1;
    tempim = double(max(tcompim(stx:edx,sty:edy,stz:edz),[],3));
    tempim = histeq(tempim/max(tempim(:)));
    iparams(:,:,:,j) = tempim;
end

%% Second selection of candidate locations of inner hair cells
Mdl = omdl2;
[pred1,scores2] = predict(Mdl,params2_2);
pred2 = double(classify(net1,iparams));     
xyzlist2 = SwitchColumn1_2(cat(1,S(idxlist).PixelList));
Isize = size(linearizedIm2);
pred2(pred1==-1)=-1;
datalist = [idxlist pred2 xyzlist scores2(:,2)];

leftend = round(min(innerCells2(:,2)));
f1 = datalist(:,4) < (leftend+200);
f2 = pred1 == 1;
leftadd = datalist(f1.*f2>0,:);

%% Detect rows of outer hair cells
disp('Estimating outer hair cells...')
for j = 1:size(datalist,1)
    tempdata = datalist(j,:);
    if tempdata(2) ~=-1
        if tempdata(4) > Isize(2)-10 || tempdata(5) < 6 || tempdata(5) > Isize(1)-6
            datalist(j,2) = -1;
        else
            tempxyz = round(tempdata(3:5));
            if tempxyz(1)+5>Isize(1)
                datalist(j,2) = -1;
            else
                
                tempint1 = linearizedIm2(tempxyz(1),tempxyz(2)+10,tempxyz(3));
                tempint2 = linearizedIm2(tempxyz(1)-5,tempxyz(2),tempxyz(3));
                tempint3 = linearizedIm2(tempxyz(1)+5,tempxyz(2),tempxyz(3));
            if tempint1 == 0 || tempint2 == 0 || tempint3 == 0
                datalist(j,2) = -1;
            end
            end
        end
    end
end

for ii = 1:3
    out1 = sortrows(datalist(datalist(:,2)==1,3:5),2);
    out2 = sortrows(datalist(datalist(:,2)==2,3:5),2);
    out3 = sortrows(datalist(datalist(:,2)==3,3:5),2);
    outs = {out1;out2;out3};
    dif_1 = diff(out1);
    dif_2 = diff(out2);
    dif_3 = diff(out3);
    difs = {dif_1;dif_2;dif_3};
    xyzlist = datalist(:,3:5);
    
    for j = 1:size(difs{ii,1})
        if difs{ii,1}(j,2) < 5
            xyz1 = outs{ii,1}(j,:);
            xyz2 = outs{ii,1}(j+1,:);
            f1 = (xyzlist(:,1) == xyz1(1)) .* (xyzlist(:,2) == xyz1(2)) .* (xyzlist(:,3) == xyz1(3));
            f2 = (xyzlist(:,1) == xyz2(1)) .* (xyzlist(:,2) == xyz2(2)) .* (xyzlist(:,3) == xyz2(3));
            if abs(xyz1(1)-xyz2(1))<4
                datalist(f1>0,2) = -1;
                datalist(f2>0,2) = -1;
            else
                datalist(f1>0,2) = 4;
                datalist(f2>0,2) = 4;
            end
        end
    end
end

f1 = datalist(:,2)==-1;
f2 = datalist(:,6)>0.8;
datalist(f1.*f2>0,2) = 4;

for ii = 1:3
    out1 = sortrows(datalist(datalist(:,2)==1,3:5),2);
    out2 = sortrows(datalist(datalist(:,2)==2,3:5),2);
    out3 = sortrows(datalist(datalist(:,2)==3,3:5),2);
    outs = {out1;out2;out3};
    dif_1 = diff(out1);
    dif_2 = diff(out2);
    dif_3 = diff(out3);
    difs = {dif_1;dif_2;dif_3};
    xyzlist = datalist(:,3:5);
    f0 = datalist(:,2) ~= -1;
    
    for j = 1:size(difs{ii,1})
        if difs{ii,1}(j,2) > 12 && difs{ii,1}(j,2) < 50
            xyz1 = outs{ii,1}(j,:);
            xyz2 = outs{ii,1}(j+1,:);
            xyz3 = interp1([xyz1(2);xyz2(2)],[xyz1;xyz2],ceil(xyz1(2)+4):floor(xyz2(2)-4));
            
            f1 = (xyzlist(:,2) > outs{ii,1}(j,2)).*(xyzlist(:,2) < outs{ii,1}(j+1,2));
            f2 = (xyzlist(:,1) > (min(outs{ii,1}(j,1),outs{ii,1}(j+1,1))-8)) ...
                .* (xyzlist(:,1) < (max(outs{ii,1}(j,1),outs{ii,1}(j+1,1))+8));
            idxs = find(f0.*f1.*f2);
            
            if ~isempty(idxs)
                idxs2 = datalist(idxs,1);
                xyz4 = SwitchColumn1_2(cat(1,S(idxs2).PixelList));
                [idxs3,d] = knnsearch(xyz4,xyz3);
                xyz5 = xyz4(idxs3(d<3,:),:);
                if ~isempty(xyz5)
                    idxs4 = sub2ind(Isize,xyz5(:,1),xyz5(:,2),xyz5(:,3));
                    idxs5 = unique(L(idxs4));
                    for kk = size(idxs5,1)
                        datalist(datalist(:,1) == idxs5(kk),2) = ii;
                        datalist(datalist(:,1) == idxs5(kk),3:5) = xyz5(1,:);
                    end
                end
            end
        end
    end
end

leftend = round(min(innerCells2(:,2)));
rightend = round(max(innerCells2(:,2)));
xyzlist = datalist(:,3:5);

interval = 60;
width = 60;
stlist = [(leftend:interval:rightend-width), rightend-width+1]';
zlist = zeros(size(stlist,1),4);
for j = 1:size(stlist,1)
    f1 = (datalist(:,2)==1).*(xyzlist(:,2) >= stlist(j)).*(xyzlist(:,2) < stlist(j)+width);
    f2 = (datalist(:,2)==2).*(xyzlist(:,2) >= stlist(j)).*(xyzlist(:,2) < stlist(j)+width);
    f3 = (datalist(:,2)==3).*(xyzlist(:,2) >= stlist(j)).*(xyzlist(:,2) < stlist(j)+width);
    temp1 = xyzlist(f1>0,3);
    temp2 = xyzlist(f2>0,3);
    temp3 = xyzlist(f3>0,3);
    zlist(j,:) = [stlist(j)+width/2 mean(temp1) mean(temp2) mean(temp3)];
end
f = sum(isnan(zlist),2)>0;
zlist(f>0,:) = [];
qzlist = interp1(zlist(:,1),zlist(:,2:4),(1:Isize(2))','linear','extrap');

for j = 1:size(datalist,1)
    if datalist(j,2)==1 || datalist(j,2)==2 || datalist(j,2)==3
        if abs(datalist(j,5)-qzlist(round(datalist(j,4)),datalist(j,2)))>9 %臒l
            datalist(j,2) = -1;
        end
    end
end

interval = 250;
width = 250;
stlist = [(leftend:interval:rightend-width), rightend-width+1]';
ydlist = zeros(size(stlist,1),2);
for j = 1:size(stlist,1)
    f1 = (datalist(:,2)==1).*(xyzlist(:,2) >= stlist(j)).*(xyzlist(:,2) < stlist(j)+width);
    f2 = (datalist(:,2)==2).*(xyzlist(:,2) >= stlist(j)).*(xyzlist(:,2) < stlist(j)+width);
    f3 = (datalist(:,2)==3).*(xyzlist(:,2) >= stlist(j)).*(xyzlist(:,2) < stlist(j)+width);
    temp1 = diff(sort(xyzlist(f1>0,2)));
    temp2 = diff(sort(xyzlist(f2>0,2)));
    temp3 = diff(sort(xyzlist(f3>0,2)));
    ydlist(j,:) = [stlist(j)+width/2 median([temp1; temp2; temp3])];
end

f = sum(isnan(ydlist),2)>0;
ydlist(f>0,:) = [];
ydlist(:,2) = medfilt1(ydlist(:,2));
qylist = interp1(ydlist(:,1),ydlist(:,2),(1:Isize(2))','linear','extrap');

for j = 1:size(datalist,1)
    if datalist(j,2)~=-1
        templabel = datalist(j,2);
        if templabel == 4
            templabel = 3;
        end
        tempy = datalist(j,4);
        templist = SwitchColumn1_2(cat(1,S(datalist(j,1)).PixelList));
        idx =  knnsearch(templist(:,3),qzlist(round(tempy),templabel));
        datalist(j,3:5) = templist(idx,:);
    end
end

for ii = 1:3
    out1 = sortrows(datalist(datalist(:,2)==1,3:5),2);
    out2 = sortrows(datalist(datalist(:,2)==2,3:5),2);
    out3 = sortrows(datalist(datalist(:,2)==3,3:5),2);
    outs = {out1;out2;out3};
    dif_1 = diff(out1);
    dif_2 = diff(out2);
    dif_3 = diff(out3);
    difs = {dif_1;dif_2;dif_3};
    xyzlist = datalist(:,3:5);
    
    for j = 1:size(difs{ii,1})
        if difs{ii,1}(j,2) < 5
            xyz1 = outs{ii,1}(j,:);
            xyz2 = outs{ii,1}(j+1,:);
            f1 = (xyzlist(:,1) == xyz1(1)) .* (xyzlist(:,2) == xyz1(2)) ...
                .* (xyzlist(:,3) == xyz1(3));
            f2 = (xyzlist(:,1) == xyz2(1)) .* (xyzlist(:,2) == xyz2(2)) ...
                .* (xyzlist(:,3) == xyz2(3));
            if abs(xyz1(1)-xyz2(1))<4
                datalist(f1>0,2) = -1;
                datalist(f2>0,2) = -1;
            else
                datalist(f1>0,2) = 4;
                datalist(f2>0,2) = 4;
            end
        end
    end
end

f1 = datalist(:,2)==-1;
f2 = datalist(:,6)>0.8;
datalist(f1.*f2>0,2) = 4;

%% Examine gaps within rows of estimated outer hair cells 
o_del = cell(3,1);
delcand = cell(3,1);
addcand = cell(3,1);
addlabel = [];
xyzlistL = GetCoordOfPositivePixels(L>0);

for ii = 1:3
    outlist = sortrows(datalist(datalist(:,2)==ii,3:5),2);
    temp = diff(outlist);
    difnorm = sqrt(sum(temp(:,1:2).^2,2));
    xyzlist = datalist(:,3:5);
    for j = 1:size(difnorm)
        if difnorm(j,1) > 11
            xyz1 = outlist(j,:);
            xyz2 = outlist(j+1,:);
            
            xyz1_3 = xyz1;
            xyz2_3 = xyz2;
            mindist = norm(xyz1-xyz2);
            
            xyz1_4 = xyz1_3;
            xyz2_4 = xyz2_3;
            dnum = max(1,round(mindist/qylist(round(xyz1(2)))-1));
            
            if dnum == 1
                xyz3 = mean([xyz1_4;xyz2_4]);
                param = ObtainPixelValuesAround2(xyz3,linearizedIm2);
                temppred = classify(net2,param);
                delcand{ii,1} = [delcand{ii,1}; xyz3];
                [param2,tempxyz,d] = ObtainFeatureQuantities(xyz3,ii,xyzlist2 ...
                    ,idxlist,L,scores2,leftend,outC);

                if temppred == '1'
                    o_del{ii,1} = [o_del{ii,1}; xyz3];
                else
                    if d < 3.5
                        tempidx = L(tempxyz(1),tempxyz(2),tempxyz(3));
                        templabel = datalist(idxlist == tempidx,2);
                        tempidx2 = find(idxlist == tempidx);
                        if templabel ~= ii
                            datalist(tempidx2,2) = ii;
                            datalist(tempidx2,3:5) = tempxyz;
                        end
                    else
                        [Lidx,Ld] = knnsearch(xyzlistL,xyz3);
                        if Ld < 3.5
                            tempxyz = xyzlistL(Lidx,:);
                            templabel = L(tempxyz(1),tempxyz(2),tempxyz(3));
                            datalist = [datalist; double(templabel),ii,tempxyz,NaN];
                            addlabel = [addlabel; [double(templabel) xyz3 score1(templabel,2)]];
                        end
                        addcand{ii,1} = [addcand{ii,1}; xyz3];
                    end
                end
            elseif dnum == 2
                vector = (xyz2_4-xyz1_4);
                vector = vector/3;
                xyz3 = [xyz1_4 + vector;xyz2_4 - vector];
                param1_1 = ObtainPixelValuesAround2(xyz3(1,:),linearizedIm2);
                param1_2 = ObtainPixelValuesAround2(xyz3(2,:),linearizedIm2);
                temppred = classify(net2,cat(4,param1_1, param1_2));
                
                delcand{ii,1} = [delcand{ii,1}; xyz3];
                
                param2_1 = ObtainFeatureQuantities(xyz3(1,:),ii,xyzlist2,idxlist,L ...
                    ,scores2,leftend,outC);
                param2_2 = ObtainFeatureQuantities(xyz3(2,:),ii,xyzlist2,idxlist,L ...
                    ,scores2,leftend,outC);
                
                for k = 1:2
                    if temppred(k,1) == '1'
                        o_del{ii,1} = [o_del{ii,1}; xyz3(k,:)];
                    else
                        [idx, d] = knnsearch(xyzlist2,xyz3(k,:));
                        if d < 3.5
                            tempxyz = xyzlist2(idx,:);
                            tempidx = L(tempxyz(1),tempxyz(2),tempxyz(3));
                            tempidx2 = find(idxlist == tempidx);
                            templabel = datalist(tempidx2,2);
                            if templabel ~= ii
                                datalist(tempidx2,2) = ii;
                                datalist(tempidx2,3:5) = tempxyz;
                            end
                        else
                            [Lidx,Ld] = knnsearch(xyzlistL,xyz3(k,:));
                            if Ld < 3.5
                                tempxyz = xyzlistL(Lidx,:);
                                templabel = L(tempxyz(1),tempxyz(2),tempxyz(3));
                                datalist = [datalist; double(templabel),ii,tempxyz, NaN];
                                addlabel = [addlabel; [double(templabel) xyz3(k,:) ...
                                    score1(templabel,2)]];
                            end
                            addcand{ii,1} = [addcand{ii,1}; xyz3(k,:)];
                        end
                    end
                end
            else 
                for kk = 1:100
                    lastxyz = xyz1_3;

                    out1 = sortrows(datalist(datalist(:,2)==1,3:5),2);
                    out2 = sortrows(datalist(datalist(:,2)==2,3:5),2);
                    out3 = sortrows(datalist(datalist(:,2)==3,3:5),2);
                    
                    [idx1] = knnsearch(out1,lastxyz,'k',2);
                    [idx2] = knnsearch(out2,lastxyz,'k',2);
                    [idx3] = knnsearch(out3,lastxyz,'k',2);
                    vector1 = out1(max(idx1),:)-out1(min(idx1),:);
                    vector1 = vector1/norm(vector1);
                    vector2 = out2(max(idx2),:)-out2(min(idx2),:);
                    vector2 = vector2/norm(vector2);
                    vector3 = out3(max(idx3),:)-out3(min(idx3),:);
                    vector3 = vector3/norm(vector3);
                    
                    vector4 = xyz2_3 - lastxyz;
                    vector4 = vector4/norm(vector4);
                    
                    vector = mean([vector1;vector2;vector3])+vector4;
                    vector(3) = 0;
                    vector = vector/norm(vector);
                    
                    nextxyz = lastxyz + vector*qylist(round(xyz1(2)));
                    d1 = norm(xyz2_3-nextxyz);
                    d2 = xyz2_3(2)-nextxyz(2);
                    if d2 < qylist(round(xyz2(2)))/2
                        break
                    elseif d1 < qylist(round(xyz2(2)))
                        nextxyz = mean([lastxyz; xyz2_3]);
                    end
                    
                    rparam2 = ObtainPixelValuesAround2(nextxyz,linearizedIm2);
                    temppred = classify(net2,rparam2);
                    
                    [rparam,tempxyz,d] = ObtainFeatureQuantities(nextxyz,ii,xyzlist2 ...
                        ,idxlist,L,scores2,leftend,outC);
                    
                    delcand{ii,1} = [delcand{ii,1}; nextxyz];
                    if temppred == '1'
                        o_del{ii,1} = [o_del{ii,1}; nextxyz];
                    else
                        if d < 3.5
                            tempidx = L(tempxyz(1),tempxyz(2),tempxyz(3));
                            templabel = datalist(idxlist == tempidx,2);
                            tempidx2 = find(idxlist == tempidx);
                            if templabel ~= ii 
                                if templabel == -1 || templabel == 4
                                    datalist(tempidx2,2) = ii;
                                    datalist(tempidx2,3:5) = tempxyz;
                                elseif ii == 1 
                                    tempxyz = lastxyz;
                                    tempxyz(1) = tempxyz(1)-5;
                                elseif ii == 2 
                                    datalist(tempidx2,2) = ii;
                                    datalist(tempidx2,3:5) = tempxyz;
                                elseif ii == 3 
                                    tempxyz = lastxyz; 
                                    tempxyz(1) = tempxyz(1)+5;
                                end
                            end
                            nextxyz = tempxyz;
                        else
                            [Lidx,Ld] = knnsearch(xyzlistL,nextxyz);
                            if Ld < 3.5
                                tempxyz = xyzlistL(Lidx,:);
                                templabel = L(tempxyz(1),tempxyz(2),tempxyz(3));
                                datalist = [datalist; double(templabel),ii,tempxyz,NaN];
                                addlabel = [addlabel; [double(templabel) nextxyz ...
                                    score1(templabel,2)]];
                                nextxyz = tempxyz;
                            end
                            addcand{ii,1} = [addcand{ii,1}; nextxyz];                          
                        end
                    end
                    xyz1_3 = nextxyz;
                end
            end
            
        end
    end
end

%% Examine the end ot rows of estimated outer hair cells 
tempadd = cell(3,1);
for j = 1:200
    out1 = sortrows([datalist(datalist(:,2)==1,3:5); tempadd{1,1}],2);
    out2 = sortrows([datalist(datalist(:,2)==2,3:5); tempadd{2,1}],2);
    out3 = sortrows([datalist(datalist(:,2)==3,3:5); tempadd{3,1}],2);
    out_lasts = [out1(end,:);out2(end,:);out3(end,:)];
    [~,lineidx] = min(out_lasts(:,2));
    lastxyz = out_lasts(lineidx,:);

    [idx1] = knnsearch(out1,lastxyz,'k',2);
    [idx2] = knnsearch(out2,lastxyz,'k',2);
    [idx3] = knnsearch(out3,lastxyz,'k',2);
    vector1 = out1(max(idx1),:)-out1(min(idx1),:);
    vector1 = vector1/norm(vector1);
    vector2 = out2(max(idx2),:)-out2(min(idx2),:);
    vector2 = vector2/norm(vector2);
    vector3 = out3(max(idx3),:)-out3(min(idx3),:);
    vector3 = vector3/norm(vector3);
    vector = mean([vector1;vector2;vector3;[0 1 0]]);
    vector = vector/norm(vector);
    
    nextxyz = lastxyz + vector*7.8;

    if nextxyz(2) > Isize(2)-10
        break
    end
    tempxyz = round(nextxyz);
    tempint = linearizedIm2(tempxyz(1),tempxyz(2)+5,tempxyz(3));
    if tempint==0
        break
    end
    
    rparam2 = ObtainPixelValuesAround2(nextxyz,linearizedIm2);
    temppred = classify(net2,rparam2);
    
    [rparam,xyz4,d] = ObtainFeatureQuantities(nextxyz,lineidx,xyzlist2,idxlist,L ...
        ,scores2,leftend,outC);
    
    delcand{lineidx,1} = [delcand{lineidx,1}; nextxyz];
    
    if temppred ~= '1' && d < 3
        nextxyz = xyz4;

        tempidx = L(xyz4(1),xyz4(2),xyz4(3));
        templabel = datalist(idxlist == tempidx,2);
        tempidx2 = find(idxlist == tempidx);
        
        if templabel ~= lineidx
            if templabel == -1 || templabel == 4
                datalist(tempidx2,2) = lineidx;
                datalist(tempidx2,3:5) = nextxyz;
            elseif lineidx == 1 
                nextxyz = lastxyz; 
                nextxyz(1) = nextxyz(1)-5;
                nextxyz(2) = nextxyz(2)+1;
                tempadd{lineidx,1} = [tempadd{lineidx,1}; nextxyz];
            elseif lineidx == 2 
                datalist(tempidx2,2) = lineidx;
                datalist(tempidx2,3:5) = nextxyz;
            elseif lineidx == 3 
                datalist(tempidx2,2) = lineidx;
                datalist(tempidx2,3:5) = nextxyz;
            end
        end
        
    elseif temppred == '1' 
        o_del{ii,1} = [o_del{ii,1}; nextxyz];
        tempadd{lineidx,1} = [tempadd{lineidx,1}; nextxyz];
    else
        [Lidx,Ld] = knnsearch(xyzlistL,nextxyz);
        if Ld < 3 
            tempxyz = xyzlistL(Lidx,:);
            templabel = L(tempxyz(1),tempxyz(2),tempxyz(3));
            datalist = [datalist; double(templabel),lineidx,tempxyz,NaN];
            addlabel = [addlabel; [double(templabel) tempxyz score1(templabel,2)]];
        else
            tempadd{lineidx,1} = [tempadd{lineidx,1}; nextxyz];
        end
    end
end

out1 = sortrows(datalist(datalist(:,2)==1,3:5),2);
out2 = sortrows(datalist(datalist(:,2)==2,3:5),2);
out3 = sortrows(datalist(datalist(:,2)==3,3:5),2);
outlist = [out1; out2; out3];
out4 = sortrows(datalist(datalist(:,2)==4,3:5),2);
delcandlist = [cat(1,delcand{:}); out4];

add4 = [];
for k = 1:20
    if isempty(add4)
        [~,d] = knnsearch(outlist(:,1:2),out4(:,1:2));
    else
        [~,d] = knnsearch([outlist(:,1:2); add4(:,1:2)],out4(:,1:2));
    end
    f = (d>6).*(d<10);
    add4 = [add4; out4(f>0,:)];
end
f = abs(add4(:,3) - qzlist(round(add4(:,2)),3))>5;
add4(f>0,:) = [];

temp = [];
for k = 1:size(leftadd,1)
    f = datalist(:,1)==leftadd(k,1);
    if datalist(f>0,2)==-1 || datalist(f>0,2)==4
        [~,d] = knnsearch([outlist(:,1:2); add4(:,1:2)],datalist(find(f),3:4));
        if d>5
           temp = [temp; find(f)];
        else
            d;
        end
    end
end

add5 = datalist(temp,3:5);
f = abs(add5(:,3) - mean(qzlist(round(add5(:,2)),:),2))>5;
add5(f>0,:) = [];
add4 = unique([add4;add5],'rows');

tempparam = zeros(69,39,1,size(add4,1));
for i = 1:size(add4,1)
    tempparam(:,:,:,i) = ObtainPixelValuesAround2(add4(i,:),linearizedIm2);
end
temppred = classify(net2,tempparam);
add4(temppred=='1',:)=[];

outlist = [sortrows(datalist(datalist(:,2)==1,3:5),2); sortrows(datalist(datalist(:,2) ...
    ==2,3:5),2); sortrows(datalist(datalist(:,2)==3,3:5),2); add4];

incim = DrawMarkedIm(linearizedIm2,round(outlist),1,1);
incim = uint16(incim*1000);
ImWrite3D(incim,[SubFolderPath '\outerHairCells.tif']);
xlswrite([SubFolderPath '\outerHairCells.xlsx'],sortrows(SwitchColumn1_2(round(outlist))),1)
save([SubFolderPath '\data.mat'],'outlist','-append')
disp('outer hair cells estimated!')
