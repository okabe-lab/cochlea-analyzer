function z = AdjustZCoord(compim,point,w,interval)
%*****************************************************************************
%This function uses "compim", linearized image of organ of corti, "point",
%specify calculation point, "w", width which specify range of calculation,
%"interval", lattice size. The output is  "z", adjusted z coordinate.
%*****************************************************************************

%% Calculate adjusted z coordinate
Isize = size(compim);
tempim = compim(:,max(1,point(2)-w):min(Isize(2),point(2)+w),:);
tempmax = max(tempim,[],3);
tempx = (interval/2:interval:size(tempmax,1)-interval/2)';
tempy = (interval/2:interval:size(tempmax,2)-interval/2)';

zs = zeros(size(tempx,1)*size(tempy,1),1);
c = 1;
for j = 1:size(tempx,1)
    for k = 1:size(tempy,1)
        tempim2 = tempim(tempx(j,1)-interval/2+1:tempx(j,1)+...
            interval/2,tempy(k,1)-interval/2+1:tempy(k,1)+interval/2,:);
        tempint = sum(sum(tempim2));
        tempint = tempint(:);
        f = find(tempint > prctile(tempint,80));
        if isempty(f)
            f = 0;
        end
        zs(c) = f(1);
        c = c+1;
    end
end
zs(zs==0)=[];
z = mean(zs);