function [ D1VarDxn1 ] = d1_dxn1_c( cmat1, var, n )
%D1_DXN1_C Compute the 1st-order derivative of variable
%   The 1st-order derivative of variable is computed using the correponding
%   coefficient matrices
%
% 	cmat1   - 1st coefficient matrices
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
   	D1VarDxn1 = fft_backward(cmat1.x1(:,1:vsize(2),1:vsize(3)).*var);
    
  case 2
    D1VarDxn1 = fft_backward(cmat1.x2(1:vsize(1),:,1:vsize(3)).*var);
    
  case 3
    D1VarDxn1 = fft_backward(reshape(...
      reshape(var,[],vsize(3))*transpose(cmat1.x3),...
      vsize(1),vsize(2),vsize(3)));
    
  otherwise
    warning('The index of direction is supported!')
  
end



end

