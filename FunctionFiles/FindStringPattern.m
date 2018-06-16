function containPattern = FindStringPattern(stringCell,stringPattern)
%*****************************************************************************
%This function uses "stringCell", cell array of strings, and "stringPattern",
%string pattern to find. The output is "conteinPattern", 1 or 0.
%*****************************************************************************

x = strfind(stringCell,stringPattern,'ForceCellOutput',true);
containPattern = not(cellfun('isempty',x));
