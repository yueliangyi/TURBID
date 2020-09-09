function [ D1VarDxn1 ] = d1_dxn1_c( cmat1, var, n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% the input data is in 3D cplx space!

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

