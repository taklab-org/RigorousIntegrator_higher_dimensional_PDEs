clear
clc


%--------------------------------------------------------------------------
%load ba_2D_DC_alt.mat
 load ba_2D_DC.mat
% Parameters
% -------------------------------------------------------------------------
Nop_coeff = [0,0,0,0;0,0,0,0;0,0,0,-1];
 h = 2^-7;
nu = 1.4;
m = [2,2];
m1m2 = prod(m);
% n = size(ba,1);
% padding for success
[n,N1,N2] = size(ba);
n_pad = 1;
ba_ext = zeros(n+n_pad,N1,N2);
ba_ext(1:n,:,:) = ba;
ba = ba_ext;
n = n + n_pad;
L = [1 1.1];
% Variation Problem
% ------------------------------------------------------------------------
[Phi,r1,Y_01,Z_01,Z_11] = rig_2D_diblock_copolymer_variational(ba,m,Nop_coeff,h,epsilon,sigma,nu,L,lambda);
Phi = reshape(Phi,n,m1m2,m1m2);
Phi(2:end,:,:) = 2*Phi(2:end,:,:);
% Adjoint Problem
nu_var = 1.65;
% ------------------------------------------------------------------------
[Psi,r2,Y_02,Z_02,Z_12] = rig_2D_diblock_copolymer_adjoint(ba,m,Nop_coeff,h,epsilon,sigma,nu,L,lambda);
Psi = reshape(Psi,n,m1m2,m1m2);
Psi(2:end,:,:) = 2*Psi(2:end,:,:);

% 
permute(sum(Phi,1),[2 3 1])*permute(sum(Psi,1),[2 3 1]) % this should be identity matrix
norm(eye(m(1)*m(2),m(1)*m(2))-permute(sum(Phi,1),[2 3 1])*permute(sum(Psi,1),[2 3 1]),1)



return
% Loading data
%--------------------------------------------------------------------------
load ba_2D_DC.mat
% Parameters
% profile on
% -------------------------------------------------------------------------
  %Nop_coeff = [0,0,0,0;0,0,0,0;0,0,0,-1];
 Nop_coeff = [0,0,0,-1];
h = 2^-7;
nu = 1.3;
%epsilon = 1;
% sigma = 0.1;
L = [1 1.1 1.2];
m = [2,2];
m1m2 = prod(m);
n = size(ba,1);

e_j = zeros(m);
e_j(4) = 1;

a_bar = permute(ba,[2,3,1]);

c = rand([m,n(1)]);
mu_k = mu_k_dcCH(size(c)-1,epsilon,sigma,L);

% a_bar = rand(5,3);
% c = rand(size(a_bar));
% e_j = zeros(5,1); e_j(1) = 1;
% mu_k = ones(size(c));
lambda = 0;
[c,~] = newton_F(mu_k,Nop_coeff,a_bar,c,h,e_j,L);
% F_lin_fin_dim(mu_k,Nop_coeff,a_bar,c,h,e_j)
% F = F_lin_fin_dim(mu_k,Nop_coeff,a_bar,c,h,e_j,L);



return
DF = DF_lin_fin_dim(mu_k,Nop_coeff,a_bar,h,c,L);
norm_1_nu_int(intval(a_bar),nu)^2./(nu.^(1:10))
% DF_fin = fin_dif(mu_k,Nop_coeff,a_bar,c,e_j,L);
% % [c,~] = newton_F(mu_k,Nop_coeff,a_bar,c,h,e_j);
% % F_lin_fin_dim(mu_k,Nop_coeff,a_bar,c,h,e_j)
% norm(reshape(DF,[],1)-reshape(DF_fin,[],1))
% 
% aa = convapbqcos(a_bar,1,a_bar,1);
% D_000_F_001 = h/2*(-mu_k(1,1,1) - 3*aa(1,1,3) + 3*aa(1,1,1));

% 
% Ma = size(a_bar);
% c_pad = zeros(Ma+1);
% c_pad(1:2,1:2,1:6) = c;
% aac = convapbqcos(a_bar,2,c_pad,1);
% % F_002 = 2*2*c(1,1,3) + h/2*( mu_k(1,1,3)*c(1,1,4) - mu_k(1,1,3)*c(1,1,2) - 3*aac(1,1,4) + 3*aac(1,1,2)  )

% 
% profile viewer
% profile off
% 







