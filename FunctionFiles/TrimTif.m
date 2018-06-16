function y = TrimTif(x)
%*****************************************************************************
%This function uses "x", character string and the output is "y", character
%string.
%*****************************************************************************

%% Remove ".tif" from character string
judge = strfind(x,'.tif');
if not(isempty(judge))
    x(judge:judge+3)=[];
    y = x;
end