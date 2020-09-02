function [ data, dsize ] = r_armadillo_cplx( fileID, gridnum, noidentifier )
%R_ARMADILLO_CPLX Read the binary output by library Armadillo
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

if fileID < 0, error('Input FileID is Invalid!'); end
if nargin < 3, noidentifier = 0; end
if ~noidentifier, fgetl(fileID); end % neglect the file identifier

%%

% convert to complex storage
gridnum(1) = gridnum(1)/2+1;

% get the size of data
dsize = str2vector(fgetl(fileID));

if ((dsize(1)~=gridnum(3)) || ...
    (dsize(2)~=gridnum(2)) || ...
    mod(dsize(3),gridnum(1))~=0)
  error('Incompatible Input Domain Size!')
end

% get all the data from file as a whole cube
dcube = zeros(dsize(1)*2,dsize(2),dsize(3));
for index = 1:dsize(3)
  dcube(:,:,index) = fread(fileID,[dsize(1)*2,dsize(2)],'double');   
end
dcube = permute(complex(dcube(1:2:end,:,:),dcube(2:2:end,:,:)),[3 2 1]);

% separate the whole cube into different cells
for index = dsize(3)/gridnum(1):-1:1
  data{1,index} = dcube((1:gridnum(1))+(index-1)*gridnum(1),:,:);
end

end

