% dbstop if error
% dbclear if error
clear
% clf
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
lambda = 3;
% 
% Fourier dimensions
dim_Fourier = 15;
N1 = dim_Fourier; N2 = dim_Fourier; N = [N1,N2];

% ratio of rectangular domain
L1 = 1; L2 = 1.1; L = [L1,L2]; % ratio of rectangle
% Lq_hat = max(L);

% truncated Fourier dimensions for variational problem
m = [3,3]; % [2,2], [3,2], or greater
mu_ast = find_mu(lambda,m,L); % maximum eigenvalue in tail
tau = 2^-4; % step size of timestep
nu = 1.0; % weight for norm of sequece
tspan = [0,tau];
num_integration = 1e4;
%
u_max = lambda + .1;
%
% Initial sequece
% a0_init = zeros(N1,N2);
% Stripe pattern: equilibria1
% which_case = 'stripe';
% a0_init(1,1) = 0.2;
% a0_init(1,2) = -0.2;
% a0_init(2,1) = 0.25; 
% a0_init(1:3,1:2) = [
%    9.9118e-02  -1.3032e-01;
%    9.6626e-01   1.6457e-02;
%   -2.4086e-02   1.8339e-02]; % 217 steps to GE
% a0_init(1:3,1:2) = [
%    1.5155e-01  -1.9379e-01;
%    9.2688e-01   3.4827e-02;
%   -2.9126e-02   2.3911e-02]; % 247 steps to GE
% a0_init(1:3,1:2) = [2.4228e-01, -3.0658e-01; 8.4031e-01, 7.5658e-02; -2.5460e-02, 2.7923e-02];% 297 steps to GE
% a0_init(1:3,1:3) = [
%    2.9134e-01  -3.9120e-01   1.2088e-02;
%    7.6444e-01   1.0229e-01  -1.2358e-02;
%   -1.4485e-02   2.7790e-02   0]; % 347 steps to GE
%
% Spot pattern: equilibria2
which_case = 'spot';
% a0_init(1,1) = -0.05;
% a0_init(2,1) = 0.005;
% a0_init(2,2) = -0.15;
% Spot pattern: equilibria2
load data_spot/a_bar_SH2D_300steps.mat a_end % 12 steps to GE with tau = 2^-x
% load data_spot/a_bar_SH2D_275steps.mat a_end % 12 steps to GE with tau = 2^-x
% 
a0_init = zeros(N1,N2);
a0_init(1:size(a_end,1),1:size(a_end,2)) = a_end;
% 
% a0_init = a_end;
% 
%
err_at_endpoint_old = 0; % initial error
if ~rigorous
  err_at_endpoint = 0;
end
%
% y_local = (zeros(1,10));
% y = (zeros(num_integration,10));
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
  [ba, delta, a_end] = getting_the_solution_timestepping(N,tspan,a0_init,lambda,L,nu,rigorous);
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
    [W_tau, W_at_t1, W_J, ba_X, kappa] = verify_evolutionop(ba,tau,m,lambda,L,nu,mu_ast);
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


    Pr = @(x) B + W_Y * (L_rho(x) .* x + C) - x;

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
%     if timestep>1
%     if any(err_at_endpoint(end)./err_at_endpoint(end-1)>1.2) 
%       tau = tau/(err_at_endpoint(end)./err_at_endpoint(end-1))^3;
%       tspan(2) = tspan(1) + tau;% next time step
%       disp('adjust timestep (smaller)')
%       continue
%     elseif all(err_at_endpoint(end)./err_at_endpoint(end-1) <= 1.01) && (tau < 2^-3)% max stepsize
%       tau = tau*1.1;
%       tspan(2) = tspan(1) + tau;% next time step
%       disp('adjust timestep (larger)')
%       continue
%     end
%     end

    %%%%%%%%%% Verify Global Existence
    if verify_GE2_2DSH(sup(err_at_endpoint(end)),a_end,nu)
      disp('Global existence is proved!')
      save(['data_',which_case,'/data_tau=',num2str(tau),'_SH2D_',num2str(timestep),'steps.mat'])
      return
    end

  end
  %% Update the initial error and time interval
  %     plot_SH_profile(a_end,[L1,L2]),pause(0.01)
  %     plot_SH_profile(a_end,[L1,L2]),view([0 90]),pause(0.01)
  %     axis([0,pi,0,pi,-u_max,u_max])
  %     pause(0.01)
  if mod(timestep,5)==0
    % save data every 5 steps
    save(['data_',which_case,'/data_tau=',num2str(tau),'_SH2D_',num2str(timestep),'steps.mat'])
  end
  W_ell_prev = W_ell; W_t1_prev = W_t1;
  W_Y_prev = W_Y; W_s_prev = W_s;
%   eps_at_t1 = min(err(end-1),err(end));
  tspan(1) = tspan(2);
  tspan(2) = tspan(2) + tau;% next time step
  timestep = timestep + 1;
  a0_init = a_end;
end