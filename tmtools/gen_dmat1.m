function [ CMat1 ] = gen_dmat1( d_size, g_num )
%UNTITLED3 Generate 1st-order differential matrices
%   This function works using the grid size of spectral domain
%
%   d_size    - domain size
% 	g_num     - grid number
%   CMat1.x1  - coefficent matrix in x1
%   CMat1.x2  - coefficent matrix in x2
%   CMat1.x3  - coefficent matrix in x3
%
%==========================================================================

g_num(1) = g_num(1)/2+1;
x3 = flipud(cheb_grid(g_num(3)))';


kx1 = 2*pi/d_size(1);
ky1 = 2*pi/d_size(2);
kxn = kx1*(0:g_num(1)-1);

for iy = g_num(2)-1:-1:0
    if (iy <= g_num(2)/2) 
        kym(iy+1) = +ky1*(iy-0);
    else
        kym(iy+1) = -ky1*(g_num(2)-iy);
    end
end

CMat1.x1 = permute(repmat(reshape(complex(0,+power(kxn,1)),1,1,[]),g_num(3),g_num(2),1),[3 2 1]);
CMat1.x2 = permute(repmat(reshape(complex(0,+power(kym,1)),1,[],1),g_num(3),1,g_num(1)),[3 2 1]);
CMat1.x3 = cf_diff(x3);

end

