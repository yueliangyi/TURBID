
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

%%
aprof   = get_loflow(cfg_file.re,cfg_file.g_num(3),0);
% plot(aprof,n_x3)

%%

csize2 = g_cplx;
csize2(1) = csize2(1)/2+1;

p00var{1,4} = zeros(csize2(1),csize2(2),csize2(3))*complex(0,0);
p00var{1,3} = zeros(csize2(1),csize2(2),csize2(3))*complex(0,0);
p00var{1,2} = zeros(csize2(1),csize2(2),csize2(3))*complex(0,0);
p00var{1,1} = zeros(csize2(1),csize2(2),csize2(3))*complex(0,0);

% add mean solution
p00var{1,1}(1,1,:) = aprof;

% enhance bottom and top
p00var{1,1}(1,1,1  ) = 0;
p00var{1,1}(1,1,end) = 0;
p00var{1,2}(1,1,1  ) = 0;
p00var{1,2}(1,1,end) = 0;
p00var{1,3}(:,:,1  ) = 0;
p00var{1,3}(:,:,end) = 0;


figure
plot(squeeze(p00var{1,1}(1,1,:)),n_x3)


p00var{1,1} = permute(p00var{1,1},[3 2 1]);
p00var{1,2} = permute(p00var{1,2},[3 2 1]);
p00var{1,3} = permute(p00var{1,3},[3 2 1]);
p00var{1,4} = permute(p00var{1,4},[3 2 1]);


p00vartmp = zeros([csize2(3)*2,csize2(2),csize2(1)]);
for index = 1:4
  p00vartmp(1:2:end,1:csize2(2)/2,1:csize2(1)) = real(p00var{1,index}(:,1:csize2(2)/2,:));
  p00vartmp(2:2:end,1:csize2(2)/2,1:csize2(1)) = imag(p00var{1,index}(:,1:csize2(2)/2,:));
  p00vartmp(1:2:end,end-csize2(2)/2:end,1:csize2(1)) = real(p00var{1,index}(:,csize2(2)/2:end,:));
  p00vartmp(2:2:end,end-csize2(2)/2:end,1:csize2(1)) = imag(p00var{1,index}(:,csize2(2)/2:end,:));
  p00var{1,index} = p00vartmp;
end



w_phase_cplx(cfg_file.w_path,0,0,p00var);





