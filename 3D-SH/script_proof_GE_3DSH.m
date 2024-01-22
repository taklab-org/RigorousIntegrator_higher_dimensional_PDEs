clear
% clc
% 
rigorous = true;
% 
addpath("defect/")
addpath("defect/Chebyshev/")
% addpath('High_Dim_Cosine_Conv/')
addpath('../Code_Gab/Code Variational and Adjoint problem/')
addpath("verify_solution/")
% 
% parameters for SH
sigma = 0.04;

% Fourier dimensions
% N1 = 5; %m = N1+1;% max wave # of Fourier w.r.t x1
% N2 = 5; %n = N2+1;% max wave # of Fourier w.r.t x2
% N3 = 5; %n = N3+1;% max wave # of Fourier w.r.t x3
N1 = 4; %m = N1+1;% max wave # of Fourier w.r.t x1
N2 = 4; %n = N2+1;% max wave # of Fourier w.r.t x2
N3 = 4; %n = N3+1;% max wave # of Fourier w.r.t x3
N = [N1,N2,N3];

% ratio of prism domain
L1 = 1; L2 = 1.1; L3 = 1.2;
L = [L1,L2,L3];
% 
m = [2,1,1]; % [2,1,1], [2,2,1], [2,2,2] or greater
mu_ast = find_mu(sigma,m,L); 
% 

tau = 2^-2; % length of time step
nu = 1; % weight for norm of sequece
tspan = [0,tau];
num_integration = 1e3;
% 
% 
% 
a0_init = zeros(N1,N2,N3);
% Stripe pattern: 
which_case = 'stripe';
% a0_init(1,1,1) = 0.1;
% a0_init(1,2,1) = -.5;
a0_init(2,1,1) = -.1;
a0_init(1,1,2) = -.005;
a0_init(1,2,1) = -.005;
% 
% 
% err_at_endpoint_old = 0;
if ~rigorous
  err_at_endpoint = 0;
end
% 
W_tau_vec = intval(zeros(num_integration,1));
W_at_t1_vec = intval(zeros(num_integration,1));
W_J_vec = intval(zeros(num_integration,1));
ba_X_vec = (zeros(num_integration,1));
vepsilon_vec = intval(zeros(num_integration,1));
delta_vec = intval(zeros(num_integration,1));
tau_vec = zeros(num_integration,1);
% 
% 
timestep = 1;
W_ell_prev = []; W_t1_prev = [];
W_Y_prev = []; W_s_prev = [];
while timestep <= num_integration
  disp(['timestep: ',num2str(timestep)])
  tau_vec(timestep) = tau;
  %%%%%%%%% getting approximate solution and defect bounds
  [ba, delta, a_end] = getting_the_solution_timestepping(N,tspan,a0_init,sigma,L,nu,rigorous);
  if rigorous
    delta_vec(timestep) = delta;
    delta = sup(delta);
  end
  disp(['delta = ',num2str(delta)])
  %
  if rigorous
    vepsilon = (compute_eps(intval(ba),intval(a0_init),nu)); % tiny errors in section
    vepsilon_vec(timestep) = vepsilon;
    %
    %%%%%%%%% verify local existence within multi-steps
%     tic
%     [err,err_at_endpoint,W_at_endpoint,W_J,Wm,W_h,ba_X,kappa] = verify_local_existence(eps_all,delta,ba,h,m,lambda,L,nu);
%     toc
    [W_tau, W_at_t1, W_J, ba_X, kappa] = verify_evolutionop(ba,tau,m,sigma,L,nu,mu_ast);
    W_tau_vec(timestep) = W_tau; W_at_t1_vec(timestep) = W_at_t1; W_J_vec(timestep) = W_J; ba_X_vec(timestep) = ba_X;
    ba_X_tmp = ba_X_vec(1:timestep);

    L_rho = @(x) 3 * (2*ba_X_tmp+x) .* x; % Special for cubic nolinearity

    % W_ell
    W_ell = intval(zeros(timestep));
    W_ell(1:end-1,1:end-1) = W_ell_prev;
    W_temp = intval(ones(1,timestep));
    for i=1:(timestep-1)
      W_temp(i) = prod(W_at_t1_vec(i:timestep-1));
    end
    W_ell(timestep,:) = W_tau*W_temp;
    % W_Y
    W_Y = intval(zeros(timestep));
    W_Y(1:end-1,1:end-1) = W_Y_prev;
    W_temp = intval(ones(1,timestep));
    for i=1:(timestep-2)
      W_temp(i) = prod(W_at_t1_vec(i+1:timestep-1));
    end
    W_Y(timestep,:) = W_tau * W_temp .* (tau_vec(1:timestep))' .* [W_J_vec(1:timestep-1);1]';
    W_Y = sup(W_Y);
    B = sup(W_ell * vepsilon_vec(1:timestep));
    C = sup(delta_vec(1:timestep));


    Pr = @(x) B + W_Y * ( L_rho(x) .* x + C) - x;

    xx = verifynlss(Pr,Pr(0));
    if all(Pr((1+eps)*sup(xx))<0,1)
      err = sup(xx);
    else
      err = NaN;
      disp('proof of local existence is failed!!')
      break
    end

    % error at the endpoint
    % W_t1
    W_t1 = intval(zeros(timestep));
    W_t1(1:end-1,1:end-1) = W_t1_prev;
    W_temp = intval(ones(1,timestep));
    for i=1:(timestep-1)
      W_temp(i) = prod(W_at_t1_vec(i:timestep-1));
    end
    W_t1(timestep,:) = W_at_t1*W_temp;
    % W_s
    W_s = intval(zeros(timestep));
    W_s(1:end-1,1:end-1) = W_s_prev;
    W_temp = intval(ones(1,timestep));
    for i=1:(timestep-1)
      W_temp(i) = prod(W_at_t1_vec(i+1:timestep));
    end
    W_s(timestep,:) =  W_temp .* (tau_vec(1:timestep))' .* (W_J_vec(1:timestep))';

    err_at_endpoint = sup(W_t1 * vepsilon_vec(1:timestep) + W_s * (L_rho(err) .* err + delta_vec(1:timestep)));
    disp(['err at the endpoint: ',num2str(err_at_endpoint(end))])

    %%%%%%%%%%

%     % adjust stepsize corresponding to increase ratio of err_at_endpoint
%     if any(err_at_endpoint./err_at_endpoint_old>1.05) && timestep~=1
%       h = h/(max(err_at_endpoint./err_at_endpoint_old))^2;
%       tspan(2) = tspan(1) + h;% next time step
%       %       m = max(1, m-1);
%       disp('adjust timestep (smaller)')
%       continue
%     elseif all(err_at_endpoint./err_at_endpoint_old <= 1.01) && (h < 2.5e-3)% max stepsize
%       h = h*1.1;
%       tspan(2) = tspan(1) + h;% next time step
%       %       m = max(2, m+1);
%       disp('adjust timestep (larger)')
%       continue
%     end
    %
    if verify_GE_3DSH(sup(err_at_endpoint(end)),a_end,nu)
      disp('Global existence is proved!')
      save(['data_',which_case,'/data_tau=',num2str(tau),'_SH3D_',num2str(timestep),'steps.mat'])
      return
    end    
  end

  %% Update the initial error and time interval
%   plot_SH_isosurfaces(a_end,L),pause(0.01)
%   SaveFig(gcf,['profiles/profile_SH_3D_timestep=',num2str(timestep)])
%   axis([0,pi,0,pi,-u_max,u_max])
%   pause(0.01)
  if mod(timestep,5)==0
    % save data every 5 steps
    save(['data_',which_case,'/data_tau=',num2str(tau),'_SH3D_',num2str(timestep),'steps.mat'])
  end
  W_ell_prev = W_ell; W_t1_prev = W_t1;
  W_Y_prev = W_Y; W_s_prev = W_s;
% 
%   err_at_endpoint_old = err_at_endpoint;
  tspan(1) = tspan(2);
  tspan(2) = tspan(2) + tau;% next time step
  timestep = timestep + 1;
  a0_init = a_end;
end