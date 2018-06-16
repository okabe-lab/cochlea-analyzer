function SSresid = CircleFitRSS(point_data, center)
%*****************************************************************************
%This function uses "point_data", xy coordinate of points, "center" center
%coordinate of circle. The output is "SSresid", residual sum of squares.
%*****************************************************************************

if isempty(point_data)
    SSresid = 10000;
    return
end
point_data = point_data(:,1:2);
temp = point_data(:,1:2)-repmat(center,size(point_data,1),1);
dist = sqrt(temp(:,1).^2+temp(:,2).^2);

yresid = dist - mean(dist);
SSresid = sum(yresid.^2);
end