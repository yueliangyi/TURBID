function [ result, strcell ] = str2vector( str, delimiter )
%STR2VECTOR Convert string into a row vector with a specified delimiter
%   The input string contains several math expression separated by the 
%   specified/default delimiter. After division by the delimiter, the math
%   expression will be calculated and a result row vector will be returned.
%
%   str         - input string
% 	delimiter   - default as space
% 	result      - result row vector
% 	strcell     - expression in cell
%
%==========================================================================

if nargin < 2 
  delimiter = ' ';
end

strcell = strsplit(str,delimiter);
for index = length(strcell):-1:1
  result(index) =  eval(strcell{index});
end

end

