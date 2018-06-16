function line = RemoveOutliersFromRow(line, num, th)
%*****************************************************************************
%This function uses "line", xyz coordinates of points, "num", number of 
%points for each calculation, "th" threshold for outlier. The output is
%"line", set of xyz coordinates of selected points.
%*****************************************************************************

%% Remove outlier to make a single row of points 
for j = 1:200
    dif1 = zeros(size(line,1),1);
    for i = 1:size(line)
        temppoint = line(i,:);
        st = max(1,min(max(1,i-floor(num/2)),size(line,1)-num+1));
        ed = min(st+num-1,size(line,1));
        temppoint2 = line(st:ed,:);
        temppoint2(temppoint2(:,2)==temppoint(2),:) = [];
        X = [ones(length(temppoint2(:,2)),1) temppoint2(:,2)];
        b = X\temppoint2(:,1);
        dif1(i,:) = abs(temppoint(1) - [1 temppoint(2)]*b);
    end
    if max(dif1)<th
        break
    else
        [~,idx2] = max(dif1);
        line(idx2,:) = [];
    end
    if size(line,1) < num
        break
    end
end