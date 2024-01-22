close all
clear
clc


% Loading data

%--------------------------------------------------------------------------
load ba_2D_DC_alt.mat
% load ba_2D_DC.mat
% Parameters
% -------------------------------------------------------------------------
Nop_coeff = [0,0,0,0;0,0,0,0;0,0,0,-1];
% h = 2^-7;
nu = 1.4;
m = [5,4];
m1m2 = prod(m);
% n = size(ba,1);
% padding for success
[n,N1,N2] = size(ba);
n_pad = 10;
ba_ext = zeros(n+n_pad,N1,N2);
ba_ext(1:n,:,:) = ba;
ba = ba_ext;
n = n + n_pad;
L = [1 1.1];
lambda = 0;
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