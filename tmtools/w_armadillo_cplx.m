function [ data, datasize ] = w_armadillo_cplx( fileID, data )
%W_ARMADILLO_CPLX Read the binary output by library Armadillo
%   Armadillo is a high quality linear algebra library (matrix maths) for
%   the C++ language, aiming towards a good balance between speed and ease
%   of use (http://arma.sourceforge.net/).
%
%   The target binary file is created by the function "save()" in Armadillo
%   with the data type of "arma_binary" which includes the size of the 
%   stored data listed as the second line.
%
%   fileID      - file identifier
%   data        - result data
%   dsize       - the size of data
%
%   The structure of data
%
%                  z
%                 ^
%                /
%               /
%              /
%            0 ---------> x
%             |
%             |  bottom
%             |
%             v
%             y
%
%==========================================================================

if fileID < 0
  error('Input FileID is Invalid!'); 
end

datasize = size(data);
datasize(1) = datasize(1)/2;

fprintf(fileID,'%s\n','ARMA_CUB_BIN_FC016');
fprintf(fileID,'%d %d %d\n',datasize);
fwrite(fileID,data,'double');

end

