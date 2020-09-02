function [ FVar ] = fft_forward( Var )
%FFT_FORWARD Fast discrete Fourier transform (DFT)
%   Transform the input variable from real space to wave number domain, but
%   in horizontal directions when the variable is 2D or 3D. The function
%   also takes care of situations with 1D inputs. In such situations, the
%   data size denotes the direction where the DFT will be take. For
%   example, a DFT in x1 direction will be returned if the input variable
%   having a size of [sz 1]. When the input has a size of [1 sz], the DFT
%   will be taken in x2 direction, instead. If the data is actually 2D but
%   the DFT in one direction is wanted, the input variable should has three
%   dimensions but keep one horizontal dimension as one.
%
%==========================================================================

DataSize    = size(Var);                    % input data size
FwdFFTScale = DataSize(1)*DataSize(2);      % scale introduced in DFT
FVar        = fft2(Var)/FwdFFTScale;        % apply DFT

% Apply truncation in x1 direction if necessary
if DataSize(1) > 1
    FVar = FVar(1:DataSize(1)/2+1,:,:);
end

end

