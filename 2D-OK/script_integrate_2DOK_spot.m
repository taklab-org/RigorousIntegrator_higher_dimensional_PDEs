dbstop if error
% dbclear if error
% 
clear
% clf

% clc
% 
rigorous = true;
% 
addpath("defect/")
addpath("defect/Chebyshev/")
addpath('High_Dim_Cosine_Conv/')
addpath('../Code_Gab/Code Variational and Adjoint problem/')
addpath("verify_solution/")

% parameters for DC
lambda = 0; % average of solution
sigma = 1; % sigma = 0 is Cahn-Hilliard
% epsilon = 1/sqrt(15); % diffusion
epsilon = 1/2.5; % diffusion
% epsilon = 1/7; % diffusion
% epsilon = 1/12; % diffusion
params = [lambda,sigma,epsilon];

% Fourier dimensions
% N1 = 24; %m = N1+1;% max wave # of Fourier w.r.t x1
% N2 = 24; %n = N2+1;% max wave # of Fourier w.r.t x2
N1 = 12; %m = N1+1;% max wave # of Fourier w.r.t x1
N2 = 12; %n = N2+1;% max wave # of Fourier w.r.t x2
% N1 = 6; %m = N1+1;% max wave # of Fourier w.r.t x1
% N2 = 6; %n = N2+1;% max wave # of Fourier w.r.t x2
N = [N1,N2];

% ratio of rectangular domain
L1 = 1; L2 = 1.1; L = [L1,L2];

% truncated Fourier dimensions for variational problem
m = [3,3];
mu_ast = find_mu(params,m,L); % maximum eigenvalue in tail
% 
tau = 2^-2; % length of time step
nu = 1.001; % weight for norm of sequece
tspan = [0,tau];
% 
% u_max =1; % works?
% 
num_integration = 1e3;
% 
% initial sequences (including the zero mode)
% a0_init = 1e-16*randn(N1,N2);
a0_init = zeros(N1,N2);
% stripe pattern: equilibria1
% a0_init(3,1) = .002;
% a0_init(1,2) = .02;
% which_case = 'stripe';
% spot pattern: equilibria2
a0_init(2,1) = -.02;
a0_init(2,2) = .001;
which_case = 'spot';
% 
a0_init(1,1) = lambda; % constraint for average of solution
% 
% err_at_endpoint_old = 0;
% 
% y_local = (zeros(1,10));
% y = (zeros(num_integration,10));
W_tau_vec = intval(zeros(num_integration,1));
W_tau_q_vec = intval(zeros(num_integration,1));
W_J_vec = intval(zeros(num_integration,1));
W_J_q_vec = intval(zeros(num_integration,1));
W_at_t1_vec = intval(zeros(num_integration,1));
% 
ba_norm_vec = (zeros(num_integration,1));
vepsilon_vec = intval(zeros(num_integration,1));
delta_vec = intval(zeros(num_integration,1));
tau_vec = zeros(num_integration,1);
% 
% 
timestep = 1;
W_ell_prev = []; W_t1_prev = [];
W_X_prev = []; W_s_q_prev = [];
W_Y_prev = []; W_s_prev = [];
while timestep <= num_integration
  disp(['timestep: ',num2str(timestep)])
  tau_vec(timestep) = tau;
  %%%%%%%%% getting approximate solution and defect bounds
  [ba, delta, a_end] = getting_the_solution_timestepping(N,tspan,a0_init,params,L,nu,rigorous);
  if rigorous
    delta_vec(timestep) = delta;
    delta = sup(delta);
  end
  disp(['delta = ',num2str(delta)])
  %
  if rigorous
    vepsilon = compute_eps(intval(ba),intval(a0_init),nu); % tiny error in section
    vepsilon_vec(timestep) = vepsilon;
    %
    %%%%%%%%% verify local existence within multi-steps
    [W_tau, W_tau_q, W_J, W_J_q, W_at_t1, ba_norm, kappa, gam] = verify_evolutionop(ba,tau,m,params,L,nu,mu_ast);
    % return
    W_tau_vec(timestep) = W_tau; W_tau_q_vec(timestep) = W_tau_q;
    W_J_vec(timestep) = W_J;         W_J_q_vec(timestep) = W_J_q;
    W_at_t1_vec(timestep) = W_at_t1; 
    ba_norm_vec(timestep) = ba_norm;
    ba_X_tmp = ba_norm_vec(1:timestep);

    L_rho = @(x) 3 * (2*ba_X_tmp+x) .* x; % Special for cubic nolinearity

    % W_ell
    W_ell = intval(zeros(timestep));
    W_ell(1:end-1,1:end-1) = W_ell_prev;
    W_temp = intval(ones(1,timestep));
    for i=1:(timestep-1)
      W_temp(i) = prod(W_at_t1_vec(i:timestep-1));
    end
    W_ell(timestep,:) = W_tau*W_temp;
    % 
    % W_X
    W_X = intval(zeros(timestep));
    W_X(1:end-1,1:end-1) = W_X_prev;
    W_temp = intval(ones(1,timestep-1));
    for i=1:(timestep-2)
      W_temp(i) = prod(W_at_t1_vec(i+1:timestep-1));
    end
    W_X(timestep,1:end-1) = W_tau * W_temp .* (tau_vec(1:timestep-1).^(1-gam)./(1-gam))' .* W_J_q_vec(1:timestep-1)';
    W_X(timestep,end) = W_tau_q * (tau_vec(timestep).^(1-gam)./(1-gam));
    W_X = sup(W_X);
    % 
    % W_Y
    W_Y = intval(zeros(timestep));
    W_Y(1:end-1,1:end-1) = W_Y_prev;
    W_temp = intval(ones(1,timestep));
    for i=1:(timestep-2)
      W_temp(i) = prod(W_at_t1_vec(i+1:timestep-1));
    end
    W_Y(timestep,:) = W_tau * W_temp .* (tau_vec(1:timestep))' .* [W_J_vec(1:timestep-1);1]';
    % W_Y = sup(W_Y);
    B = sup(W_ell * vepsilon_vec(1:timestep));
    C = sup(W_Y * delta_vec(1:timestep));

    Pr = @(x) B + W_X * (L_rho(x) .* x) + C - x;

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
    % 
    % W_s_q
    W_s_q = intval(zeros(timestep));
    W_s_q(1:end-1,1:end-1) = W_s_q_prev;
    W_temp = intval(ones(1,timestep));
    for i=1:(timestep-1)
      W_temp(i) = prod(W_at_t1_vec(i+1:timestep));
    end
    W_s_q(timestep,:) =  W_temp .* (tau_vec(1:timestep).^(1-gam)./(1-gam))' .* (W_J_q_vec(1:timestep))';
    % 
    % W_s
    W_s = intval(zeros(timestep));
    W_s(1:end-1,1:end-1) = W_s_prev;
    W_temp = intval(ones(1,timestep));
    for i=1:(timestep-1)
      W_temp(i) = prod(W_at_t1_vec(i+1:timestep));
    end
    W_s(timestep,:) =  W_temp .* (tau_vec(1:timestep))' .* (W_J_vec(1:timestep))';

    err_at_endpoint = sup(W_t1 * vepsilon_vec(1:timestep) + W_s_q * (L_rho(err) .* err) + W_s * delta_vec(1:timestep));
    disp(['err at the endpoint: ',num2str(err_at_endpoint(end))])

    %%%%%%%%%%

    % adjust stepsize corresponding to increase ratio of err_at_endpoint
    if timestep>1
    if any(err_at_endpoint(end)./err_at_endpoint(end-1)>1.2) 
      tau = tau/(err_at_endpoint(end)./err_at_endpoint(end-1))^3;
      tspan(2) = tspan(1) + tau;% next time step
      disp('adjust timestep (smaller)')
      continue
    elseif all(err_at_endpoint(end)./err_at_endpoint(end-1) <= 1.01) && (tau < 2^-2)% max stepsize
      tau = tau*1.1;
      tspan(2) = tspan(1) + tau;% next time step
      disp('adjust timestep (larger)')
      continue
    end
    end
  % 
  if mod(timestep,5)==0
    % save data every 5 steps
    save(['data_',which_case,'/data_tau=',num2str(tau),'_OK2D_',num2str(timestep),'steps.mat'])
  end
  W_ell_prev = W_ell; W_t1_prev = W_t1;
  W_X_prev = W_X; W_s_q_prev = W_s_q;
  W_Y_prev = W_Y; W_s_prev = W_s;
  end

  %% Update the initial error and time interval
  plot_DC_profile(a_end,L),view([0 90]),pause(0.01)
%   plot_DC_pattern(a_end,L),view([0 90]),pause(0.01)
%   plot_DC_profile(a_end,L),pause(0.01)
%   axis([0,pi,0,pi,-u_max,u_max])
  tspan(1) = tspan(2);
  tspan(2) = tspan(2) + tau;% next time step
  timestep = timestep + 1;
  a0_init = a_end;
end