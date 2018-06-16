function list = DetectSingleRowOfPoints(current, candidate, xdist, mindist)
%*****************************************************************************
%This function uses "current", xy coordinates of start point for search,
%"candidate", xy coordinates of points, "xdist", maximum distance in x
%direction for searching neighbor points, and "mindist", minimum distance 
%in y direction for searching neighbor points. The output is "list", xy 
%coordinates of selected points.
%*****************************************************************************

%% Extract a single row of points containing specified point
    function [newlist,candid] = search(point, candid, list, tf)
        xdist_list = abs(candid(:,1)-point(1));
        f1 = (xdist_list < xdist).*(abs(candid(:,2)-point(2))<32);
        if tf == 1
            f2 = candid(:,2) > (point(2) + mindist);
        else
            f2 = candid(:,2) < (point(2) - mindist);
        end
        if sum(f2)==0
            newlist = list;
            return
        end
        
        candid2 = candid;
        candid2(f1.*f2==0,:)=[];
        if isempty(candid2)
            candid2 = candid(f2>0,:);
        end
        dist = sqrt((candid2(:,1)-point(1)).^2 + (candid2(:,2)-point(2)).^2);
        [~,index] = min(dist);
        point = candid2(index(1),:);
        
        f3 = abs(candid2(:,2) - point(2))< 6;
        point2 = candid2(f3>0,:);
        if sum(abs(point2(:,1)-point(1))<=2)>0
            [~,idx] = min(abs(point2(:,1)-point(1)));
        else
            [~,idx] = max(point2(:,1));
        end
        point = point2(idx,:); 

        list = [list; point];
        [newlist, candid] = search(point,candid,list,tf);
        
        return
    end

index2 = find((candidate(:,1) == current(1)).*(candidate(:,2) == current(2)));
candidate(index2(1),:) = [];

list1 =  current;
[list1,~] = search(current, candidate,list1,1);

list2 = current;
[list2,~] = search(current, candidate,list2,0);

list = [flip(list2(2:end,:,:),1); list1];
end