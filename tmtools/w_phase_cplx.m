function W_PHASE_CPLX( workpath, folderindex, phaseindex, p0Xvar )
%READPHASEVARIABLE Read the binary output of phase variable
%   Written in binary, the target file may contain output of several
%   variables with different IO level. The corresponding sequence of these
%   variables in the target file can be found in program code.
%
%   location    - folder contains files
%   indstep     - index of step
%   indphase    - index of phase
%   extend      - extend the x1-x2 plane (default no)
%   data        - result cell of data
%   dsize       - typical size of data member
%
%==========================================================================

if nargin < 3
  error('Not Enough Input Arguments!'); 
end

location = fullfile(workpath,num2str(folderindex));
if ~exist(location,'dir')
  mkdir(location);
end


filename = ['phase_',sprintf('%02d',phaseindex),'_variable.dat'];
filepath = fullfile(location,filename);
fileID = fopen(filepath,'w+');
if (fileID <= 0)
  error(['Can not Create File ' filepath '!']); 
end


% add the label
fprintf(fileID,'%s\n','IOControlLevelOne');


if phaseindex == 0
  w_armadillo_cplx(fileID,cat(3,p0Xvar{1,1},p0Xvar{1,2},p0Xvar{1,3}));
  w_armadillo_cplx(fileID,p0Xvar{1,4});
else
  w_armadillo_cplx(fileID,p0Xvar{1,1});
end

fclose(fileID);


end

