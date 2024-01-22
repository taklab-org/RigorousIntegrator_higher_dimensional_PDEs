function [W_h, W_h_at_t1, W_J, ba_norm, kappa] = verify_evolutionop(ba,h,m,lambda,L,nu,mu_ast)
% [err,err_at_endpoint,,W_J,Wm,W_h,ba_X,kappa]  = 
% 
[n,N1,N2,N3] = size(ba);
% n = size(ba,1);
%% estimate the evolution operator
ba_norm = wnorm(reshape(chebmag(reshape(intval(ba),n,N1*N2*N3)),N1,N2,N3),nu);
% ba_X_alt = wnorm(permute(sum(abs(intval(ba)),1),[2,3,4,1]),nu);
%
[~,ba2]=powerconvcos(intval(ba),2);
size_ba2 = size(ba2);
n_ba2 = size_ba2(1); N1_ba2 = size_ba2(2); N2_ba2 = size_ba2(3); N3_ba2 = size_ba2(4);
ba2_sup = reshape(chebmag(reshape(ba2,n_ba2,N1_ba2*N2_ba2*N3_ba2)),N1_ba2,N2_ba2,N3_ba2);
% ba2_sup_alt = permute(sum(abs(ba2),1),[2,3,4,1]);
ba2_norm = wnorm(ba2_sup,nu);
%
% mu_ast = intval(lambda) - (1-(m(1)*intval(L(1))).^2-(m(2)*intval(L(2))).^2-(m(3)*intval(L(3))).^2).^2;
%
%
% setting for SH
g_ba = 3*ba2_norm; % g(a_bar)
% Lq_hat = max(L);
cm = getting_cm(ba2_sup,m,nu);
cinf = getting_cinf(ba2_sup,m,nu);
ba2_sup(1,1)=0; c_conv=3*wnorm(ba2_sup,nu);
cm = min(cm,c_conv);
cinf = min(cinf,c_conv);
vtheta = mu_ast+g_ba;
%
% Bounds for the evolution operator in the tail
W_infinity = (exp(vtheta*h)-1)/vtheta;
W_bar_infty = (W_infinity-h)/vtheta;
W_infinite_sup = max(1,exp(vtheta*h));

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
    exp(vtheta*h) +   Wm*W_bar_infty*cm*cinf*U_matrix(2,2)];
  U_matrix_J = [Wm_J+Wm_J*W_bar_infty*cm*cinf*U_matrix(1,1),...
    Wm_J*cm*W_infinity +  Wm_J*W_bar_infty*cm*cinf*U_matrix(1,2);...
    Wm*W_infinity*cinf + Wm*W_bar_infty*cm*cinf*U_matrix(2,1),...
    W_infinite_sup +   Wm*W_bar_infty*cm*cinf*U_matrix(2,2)];
  W_h = (norm(U_matrix,1)); % Matrix 1-norm
  ba_norm = sup(ba_norm);% ba2_X = sup(ba2_X);
  W_h_at_t1 = min(W_h,(norm(U_matrix_at_t1,1)));
  W_J = min(W_h,sup(norm(U_matrix_J,1)));
else
  disp('Linearized problem is not solved (kappa<0)')
  W_h=NaN; W_h_at_t1=NaN; W_J=NaN;
  return
end
%