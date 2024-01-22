close all
clear
clc


% Loading data

%--------------------------------------------------------------------------
load ba_3D_SH.mat
 profile on
% Parameters
% -------------------------------------------------------------------------
tic
ba = ba(:,1:3,1:3,1:3);
Nop_coeff = [0,0,0,-1];
h = 2^-7;
nu = 1.3;
lambda = 3;
m = [2 2 2];
m1m2m3 = prod(m);
n = size(ba,1);
L = [1 1 1 1];
% Variation Problem
% ------------------------------------------------------------------------
[Phi,r1,Y_01,Z_01,Z_11] = rig_3D_SH_variational(ba,m,Nop_coeff,h,lambda,nu,L);
Phi = reshape(Phi,n,m1m2m3,m1m2m3);
Phi(2:end,:,:,:) = 2*Phi(2:end,:,:,:);
% Adjoint Problem
% ------------------------------------------------------------------------
[Psi,r2,Y_02,Z_02,Z_12] = rig_3D_SH_adjoint(ba,m,Nop_coeff,h,lambda,nu,L);
Psi = reshape(Psi,n,m1m2m3,m1m2m3);
Psi(2:end,:,:,:) = 2*Psi(2:end,:,:,:);

% 
permute(sum(Phi,1),[2 3 4 1])*permute(sum(Psi,1),[2 3 4 1]) % this should be identity matrix
norm(eye(m(1)*m(2)*m(3),m(1)*m(2)*m(3))-permute(sum(Phi,1),[2 3 4 1])*permute(sum(Psi,1),[2 3 4 1]),1)
toc
profile viewer
profile off