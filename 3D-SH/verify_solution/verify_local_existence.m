function [err,err_at_endpoint,W_at_endpoint,W_J,Wm,W_h,ba_X,kappa] = verify_local_existence(eps_all,delta,ba,h,m,lambda,L,nu)
% 
[n,N1,N2,N3] = size(ba);
%% estimate the evolution operator
ba_X = wnorm(reshape(chebmag(reshape(intval(ba),n,N1*N2*N3)),N1,N2,N3),nu);
% ba_X_alt = wnorm(permute(sum(abs(intval(ba)),1),[2,3,4,1]),nu);
%
[~,ba2]=powerconvcos(intval(ba),2);
size_ba2 = size(ba2);
n_ba2 = size_ba2(1); N1_ba2 = size_ba2(2); N2_ba2 = size_ba2(3); N3_ba2 = size_ba2(4);
ba2_sup = reshape(chebmag(reshape(ba2,n_ba2,N1_ba2*N2_ba2*N3_ba2)),N1_ba2,N2_ba2,N3_ba2);
% ba2_sup_alt = permute(sum(abs(ba2),1),[2,3,4,1]);
ba2_X = wnorm(ba2_sup,nu);
% ba2_X_alt = wnorm(ba2_sup_alt,nu);
%
mu_m = intval(lambda) - (1-(m(1)*intval(L(1))).^2-(m(2)*intval(L(2))).^2-(m(3)*intval(L(3))).^2).^2;
%
%
% setting for SH
g_ba = 3*ba2_X; % g(a_bar)
Lq_hat = max(L);
cm = Lq_hat * getting_cm(ba2_sup,m,nu);
cinf = Lq_hat * getting_cinf(ba2_sup,m,nu);
ba2_sup(1,1)=0; c_conv=3*Lq_hat*wnorm(ba2_sup,nu);
cm = min(cm,c_conv);
cinf = min(cinf,c_conv);
gamma = 1;
vtheta = mu_m+gamma*Lq_hat*g_ba;
%
W_infinity = gamma*(exp(vtheta*h)-1)/vtheta;
W_bar_infty = (W_infinity-gamma*h)/vtheta;
W_infinite_sup = gamma*max(1,exp(vtheta*h));

%% solve variational problem
Nop_coeff = [0,0,0,-1];
nu_var = 1.4;
% n = size(ba,1);
m1m2m3 = prod(m);
% 
% % padding chebyshev for variational problems
% n_pad = 0;
% ba_ext = zeros(n+n_pad,N1,N2,N3);
% ba_ext(1:n,:,:,:) = ba;
% ba = ba_ext;
% n = n + n_pad;
%
[Phi,~,Y0,Z0,Z1] = rig_3D_SH_variational(ba,m,Nop_coeff,h,lambda,nu_var,L);
r1 = Y0./(1-intval(Z0)-intval(Z1));
Phi = reshape(Phi,n,m1m2m3,m1m2m3);
%
% nu_var = 1.65;
[Psi,~,Y0,Z0,Z1] = rig_3D_SH_adjoint(ba,m,Nop_coeff,h,lambda,nu_var,L);
r2 = Y0./(1-intval(Z0)-intval(Z1));
Psi = reshape(Psi,n,m1m2m3,m1m2m3);
%
% A technique to estimate a uniform bound using interval arithmetic
[Wm,Wm_J,Wm_at_t1] = getting_Wm(Phi,Psi,r1,r2,h,n,m,nu);
Wm = sup(Wm); Wm_J = sup(Wm_J);
Phi(2:end,:,:) = 2*Phi(2:end,:,:); % convert to two-sided chebyshev
Wm_at_t1 = min(wopnorm(reshape(sum(intval(Phi),1), m1m2m3, m1m2m3)+r1, nu), Wm_at_t1);
% 
kappa = 1 - Wm*W_bar_infty*cm*cinf;
kappa = inf(kappa);
%
if kappa>0
  U_matrix = [Wm, Wm*W_infinity*cm;...
    Wm*W_infinity*cinf, W_infinite_sup]/kappa;
  U_matrix_at_t1 = [Wm_at_t1+Wm_J*W_bar_infty*cm*cinf*U_matrix(1,1),...
    Wm_J*cm*W_infinity +  Wm_J*W_bar_infty*cm*cinf*U_matrix(1,2);...
    Wm*W_infinity*cinf + Wm*W_bar_infty*cm*cinf*U_matrix(2,1),...
    gamma*exp(vtheta*h) +   Wm*W_bar_infty*cm*cinf*U_matrix(2,2)];
  U_matrix_J = [Wm_J+Wm_J*W_bar_infty*cm*cinf*U_matrix(1,1),...
    Wm_J*cm*W_infinity +  Wm_J*W_bar_infty*cm*cinf*U_matrix(1,2);...
    Wm*W_infinity*cinf + Wm*W_bar_infty*cm*cinf*U_matrix(2,1),...
    W_infinite_sup +   Wm*W_bar_infty*cm*cinf*U_matrix(2,2)];
  W_h = sup(norm(U_matrix,1)); % Matrix 1-norm
  ba_X = sup(ba_X); %ba2_X = sup(ba2_X);
  L_rho = @(x) 3 * (2*ba_X+x) .* x; % Special for cubic nolinearity
  F = @(x) W_h*(eps_all+h*(Lq_hat * L_rho(x) .* x + delta))-x;
  W_at_endpoint = min(W_h,sup(norm(U_matrix_at_t1,1)));
  W_J = min(W_h,sup(norm(U_matrix_J,1)));
else
  disp('Linearized problem is not solved (kappa<0)')
  W_at_endpoint = NaN; W_h = NaN; W_J = NaN; err = NaN;
  return
end
%
%
%% verify the contraction mapping:
[xx,xx_cand,data] = verifynlssall(F,infsup(0,1));
if ~any(xx>0)
  disp('contraction mapping is not verified!')
  err = NaN;
  return
end
while(1)
  if isempty(xx_cand)
    err = sup(xx(:,all(F((1+eps)*sup(xx))<0,1)));
    break
  else
    [xx,xx_cand,data] = verifynlssall(data);
  end
end
%% error bound at the endpoint
err_at_endpoint = intval(W_at_endpoint)*eps_all+intval(W_J)*h*(Lq_hat*L_rho(err)*err+intval(delta));
err_at_endpoint = min(err,sup(err_at_endpoint));
disp(['error at endpoint = ',num2str(err_at_endpoint)])
