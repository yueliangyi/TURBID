function [ z_l ] = cheb_grid( N_z )
%CHEB_GRID Generate chebyshev grid with a given dimension
%   The grid is of Chebyshev-Gauss-Lobatto type locating at
%                     z_l = cos(l*pi/N_z)
%   where l is the index, N_z is the total number of grid points. The
%   returned location has an interval of [-1,1] where the entry z_l(1)=-1 
%   and z_l(end)=+1
% 
%   N_z     - dimension of grid
% 	z_l     - location of gird
%
%==========================================================================

if rem(N_z,2) ~= 1
  error('Input dimension should be odd!'); 
end

z_l = cos(pi.*linspace(0,1,N_z))';

end

