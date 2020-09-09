function [ D2VarDxn2 ] = d2_dxn2_c( obj, var, n )
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
   	D2VarDxn2 = fft_backward(obj.Diff.SOrder.x1(:,1:vsize(2),1:vsize(3)).*var);
    
  case 2
    D2VarDxn2 = fft_backward(obj.Diff.SOrder.x2(1:vsize(1),:,1:vsize(3)).*var);
    
  case 3
    D2VarDxn2 = fft_backward(reshape(...
      reshape(var,[],vsize(3))*transpose(obj.Diff.SOrder.x3),...
      vsize(1),vsize(2),vsize(3)));
    
  otherwise
    warning('The index of direction is supported!')
  
end



end

