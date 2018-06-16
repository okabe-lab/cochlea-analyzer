function y = SwitchColumn1_2(x)
%*****************************************************************************
%This function uses "x", 2D array with three columns. The output is "y", 
%modified array. 
%*****************************************************************************

%% Switch column 1 and 2
    y = zeros(size(x,1),size(x,2));
    y(:,2) = x(:,1);
    y(:,1) = x(:,2);
    y(:,3) = x(:,3);