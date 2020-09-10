function [ D2VarDxn2 ] = d2_dxn2_c( cmat2, var, n )
%D2_DXN2_C Compute the 2nd-order derivative of variable
%   The 2nd-order derivative of variable is computed using the correponding
%   coefficient matrices
%
% 	cmat2   - 2nd coefficient matrices
%   var     - 3D variable in complex
%   n       - index of direction
%             1: streamwise
%             2: spanwise
%             3: vertical
%
%==========================================================================

vsize  = size(var);
if numel(size(var)) == 2
  vsize(3) = 1;
end

%%

switch n
  
  case 1
   	D2VarDxn2 = fft_backward(cmat2.x1(:,1:vsize(2),1:vsize(3)).*var);
    
  case 2
    D2VarDxn2 = fft_backward(cmat2.x2(1:vsize(1),:,1:vsize(3)).*var);
    
  case 3
    D2VarDxn2 = fft_backward(reshape(...
      reshape(var,[],vsize(3))*transpose(cmat2.x3),...
      vsize(1),vsize(2),vsize(3)));
    
  otherwise
    warning('The index of direction is supported!')
  
end



end

