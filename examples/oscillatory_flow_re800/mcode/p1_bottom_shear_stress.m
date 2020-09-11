
%
% Function: example
%   Compute the bottom shear stress for an individual output
%

clc; clear; close all

%% load the configuration file and compute the basics
cfg_file.w_path = '../result';            % working path
cfg_file.d_size = [2 2 2];          % domain size
cfg_file.g_num  = [288 1152 161];   % grid number
cfg_file.re     = 24000;            % Reynolds number

% grid points in smaller domain
cfg_file.g_cplx = [cfg_file.g_num(1:2)*2/3 cfg_file.g_num(3)];

% coefficient matrices for the 1st-order derivatives 
cmat1 = gen_dmat1(cfg_file.d_size,cfg_file.g_cplx);

% 1d and 3d grid locations
gx1_1d = linspace(0,cfg_file.d_size(1),cfg_file.g_cplx(1)+1);
gx2_1d = linspace(0,cfg_file.d_size(2),cfg_file.g_cplx(2)+1);
gx1_1d = gx1_1d(1:end-1);
gx2_1d = gx2_1d(1:end-1);
gx3_1d = flipud(cheb_grid(cfg_file.g_cplx(3)))';
[gx1_3d,gx2_3d,gx3_3d] = meshgrid(gx1_1d,gx2_1d,gx3_1d);
gx1_3d = permute(gx1_3d,[2 1 3]);
gx2_3d = permute(gx2_3d,[2 1 3]);
gx3_3d = permute(gx3_3d,[2 1 3]);

%% read the output results

% target file index
tind = 74;

% open the target file
tfloc = fullfile(cfg_file.w_path,num2str(tind));
fname = 'phase_00_variable.dat';
fpath = fullfile(tfloc,fname);
fileID = fopen(fpath);
if fileID<=0, error(['Can Not Find File ' fpath '!']); end
disp(['Read File: ' fpath]);

% neglect the marker of IO level!
fgetl(fileID);

% read velocities and pressure
p00var = r_armadillo_cplx(fileID,cfg_file.g_cplx);
p00var{1,4} = r_armadillo_cplx(fileID,cfg_file.g_cplx);
p00var{1,4} = p00var{1,4}{1,1};

% close the target file
fclose(fileID);


%%

du1dx3 = d1_dxn1_c(cmat1,p00var{1,1},3);
du2dx3 = d1_dxn1_c(cmat1,p00var{1,2},3);

tau_b = sqrt(du1dx3(:,:,1).^2+du2dx3(:,:,1).^2);


%%
close all
figure
surf(squeeze(gx1_3d(:,:,1)),squeeze(gx2_3d(:,:,1)),tau_b,'LineStyle','none')







