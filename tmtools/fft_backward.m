function [ BVar ] = fft_backward( Var )
%FFT_BACKWARD inverse discrete fast Fourier transform (IDFT)
%   Detailed explanation follows the function fft_forward() where the IDFT
%   will be taken instead, and the input variable will be transformed back
%   to the real space.
%
%==========================================================================

DataSize = size(Var);   % input data size

% Find the scale in IDFT
if DataSize(1) > 1
  BwdFFTScale = (DataSize(1)-1)*2*DataSize(2);
else
  BwdFFTScale = DataSize(2);
end

BVar = ifft2(cat(1,Var,conj(Var(end-1:-1:2,:,:))),'symmetric')*BwdFFTScale;

end

