
clc; clear; close all

%% load the configuration file and compute the basics
cfg_file.w_path = '../result';    % working path
cfg_file.g_num  = [3 3 129];      % grid number
cfg_file.re     = 180;            % Reynolds number (u_tau*H/2/nu)
cfg_file.force  = 1.0;            % dring force

% grid points in smaller domain
g_cplx = [cfg_file.g_num(1:2)*2/3 cfg_file.g_num(3)];

% find the Chebyshev-Gauss-Lobatto grid points
n_x3 = flipud(cheb_grid(cfg_file.g_num(3)))';


%% read the carrier phase and calculate the error

step    = 100;
stime   = step;
etime   = 100000;
dt      = 2*pi/12000;

nt      = etime/step;
t       = (1:nt)*dt*step;

for index = nt:-1:1
  
  % target file index
  tindex = index*step;

  % open the target file
  tfloc = fullfile(cfg_file.w_path,num2str(tindex));
  fname = 'phase_00_variable.dat';
  fpath = fullfile(tfloc,fname);
  fileID = fopen(fpath);
  if fileID<=0, error(['Can Not Find File ' fpath '!']); end
  disp(['Read File: ' fpath]);

  % neglect the marker of IO level!
  fgetl(fileID);

  % read velocities and pressure
  p00var = r_armadillo_cplx(fileID,g_cplx);
  p00var{1,4} = r_armadillo_cplx(fileID,g_cplx);
  p00var{1,4} = p00var{1,4}{1,1};

  % close the target file
  fclose(fileID);
  
  % do plane-averaging for velocities
  p00varnum = numel(p00var);
  for ind = p00varnum:-1:1
    velp(:,:,:,ind) = fft_backward(p00var{ind});
    avgs12_vp(:,ind) = squeeze(mean(mean(velp(:,:,:,ind),1),2));
  end
  
  
  a_solu = get_loflow(cfg_file.re,cfg_file.g_num(3),t(index));
  mu(:,index) = avgs12_vp(:,1);
  err(index) = sqrt(mean((a_solu-mu(:,index)).^2));
  
  
end



%% plot comparsion

close all

figure('Position',[200,300,700,400])
semilogy((1:length(err))*step*dt,err,'LineWidth',2.0)
grid on
set(gca,'LineWidth',1.0,'FontSize',12,'XTick',(0:5:50))
xlabel('dimensionless time','FontSize',14)
ylabel('root mean square error','FontSize',14)
xlim([0 50])
% export_fig ./figure/error.pdf -pdf -transparent -r1200


%%
scale = 10;
figure('Position',[200,300,700,400])
hold on
for index = scale*12:-scale:scale
    plot(mu(:,index),n_x3)
end
hold off
box on
set(gca,'LineWidth',1.0,'FontSize',12,'YTick',-1:0.2:1)
xlabel('dimensionless streamwise velocity','FontSize',14)
ylabel('dimensionless height','FontSize',14)
% export_fig ./figure/evolution.pdf -pdf -transparent -r1200










