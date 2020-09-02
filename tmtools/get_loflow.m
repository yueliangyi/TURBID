function [ a_sol ] = get_loflow( Re, N_z, t )
%GET_LOFLOW Get velocity profile of laminar oscillatory flow
%
%   Reference:
%   Duarte, A. S. R., Miranda, A. I., & Oliveira, P. J. (2008).
%   Numerical and analytical modeling of unsteady viscoelastic flows: The
%   start-up and pulsating test case problems. Journal of non-newtonian
%   fluid mechanics, 154(2-3), 153-169.
%
%   Re      - stokes Reynolds number
%   N_z     - dimension of grid
% 	t       - time in radians
%   a_sol   - analytical velocity profile
%
%==========================================================================

loc_z   = cheb_grid(N_z);
alpha   = sqrt(Re);
a_sol   = real(1i*(cosh((1+1i)/sqrt(2)*alpha*loc_z)./cosh((1+1i)/sqrt(2)*alpha)-1)*exp(1i*t));

end

