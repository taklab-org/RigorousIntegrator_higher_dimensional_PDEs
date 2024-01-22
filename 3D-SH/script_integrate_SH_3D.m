clear
% clc
% 
% success = 0;
rigorous = false;
% 
addpath("defect/")
addpath("defect/Chebyshev/")
addpath('High_Dim_Cosine_Conv/')
addpath('../Code_Gab/Code Variational and Adjoint problem/')
addpath("verify_solution/")
% 
% parameters for SH
% lambda = 7.5;
lambda = 0.04;

% Fourier dimensions
% N1 = 5; %m = N1+1;% max wave # of Fourier w.r.t x1
% N2 = 5; %n = N2+1;% max wave # of Fourier w.r.t x2
% N3 = 5; %n = N3+1;% max wave # of Fourier w.r.t x3
N1 = 6; %m = N1+1;% max wave # of Fourier w.r.t x1
N2 = 6; %n = N2+1;% max wave # of Fourier w.r.t x2
N3 = 6; %n = N3+1;% max wave # of Fourier w.r.t x3
N = [N1,N2,N3];

% ratio of prism domain
L1 = 1; L2 = 1.1; L3 = 1.2;
L = [L1,L2,L3];
% 
% m = find_m(lambda,N,L);
m = [1,1,2];
% 

h = 2^-2; % length of time step
nu = 1; % weight for norm of sequece
tspan = [0,h];
num_integration = 287;
% 
% u_max = lambda + .1;
% 
% 
a0_init = zeros(N1,N2,N3);
% a0_init(1,1,1) = 0.1;
% a0_init(1,2,1) = -.5;
a0_init(2,1,1) = -.005;
a0_init(1,1,2) = -.005;
a0_init(1,2,1) = -.005; % case 3
% 
% 
% 
err_at_endpoint_old = 0;
if ~rigorous
  err_at_endpoint = 0;
end
% 
% y_local = (zeros(1,10));
% y = (zeros(num_integration,10));
% y = []; % Data container
% 
timestep = 1;
while timestep <= num_integration
  disp(['timestep: ',num2str(timestep)])
  [ba, delta, a_end] = getting_the_solution_timestepping(N,tspan,a0_init,lambda,L,nu,rigorous);
  % if rigorous
  %   delta = sup(delta);
  % end
  disp(['delta = ',num2str(delta)])
  %
%   if rigorous
%     %   eps_all = err_at_endpoint_old;
%     eps_all = sup(err_at_endpoint_old + compute_eps(intval(ba),intval(a0_init),nu));
%     %
%     %%%%%%%%%
%     tic
%     [err,err_at_endpoint,W_at_endpoint,W_J,Wm,W_h,ba_X,kappa] = verify_local_existence(eps_all,delta,ba,h,m,lambda,L,nu);
%     toc
%     if any(isnan(err))
%       disp('local existence is failed!!')
%       break
%     end
%     %%%%%%%%%%
% 
% %     % adjust stepsize corresponding to increase ratio of err_at_endpoint
% %     if any(err_at_endpoint./err_at_endpoint_old>1.05) && timestep~=1
% %       h = h/(max(err_at_endpoint./err_at_endpoint_old))^2;
% %       tspan(2) = tspan(1) + h;% next time step
% %       %       m = max(1, m-1);
% %       disp('adjust timestep (smaller)')
% %       continue
% %     elseif all(err_at_endpoint./err_at_endpoint_old <= 1.01) && (h < 2.5e-3)% max stepsize
% %       h = h*1.1;
% %       tspan(2) = tspan(1) + h;% next time step
% %       %       m = max(2, m+1);
% %       disp('adjust timestep (larger)')
% %       continue
% %     end
% 
%     %   %% Data
%     y_local(1) = tspan(1);
%     y_local(2) = tspan(2);
%     y_local(3) = Wm;
%     y_local(4) = W_at_endpoint;
%     y_local(5) = W_h;
%     y_local(6) = ba_X;
%     y_local(7) = kappa;
%     y_local(8) = err_at_endpoint;
%     y_local(9) = delta;
%     y_local(10) = err;
%     %   y = [y;y_local];
%     y(timestep,:) = y_local;
%     %
%   end

  %% Update the initial error and time interval
  plot_SH_isosurfaces(a_end,L),pause(0.01)
  % SaveFig(gcf,['profiles/profile_SH_3D_timestep=',num2str(timestep)])

  err_at_endpoint_old = err_at_endpoint;
  tspan(1) = tspan(2);
  tspan(2) = tspan(2) + h;% next time step
  timestep = timestep + 1;
  a0_init = a_end;
end
% if rigorous
%   zero_index = find(y(2:end,1)==0,1,'first');
%   y(zero_index:end,:) = [];
% %   save(['data_h=',num2str(h),'_',which_case,'.mat'],'y')
% end
