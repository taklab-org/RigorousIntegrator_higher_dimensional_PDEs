% -------------------------------------------------------------------------
clear
clc
close all

% Loading data
%--------------------------------------------------------------------------
load ba_2D_SH.mat
% Parameters
% -------------------------------------------------------------------------
Nop_coeff = [0,0,0,-1];
h = 2^-7;
lambda = 3;
a_bar = permute(ba,[2,3,1]);
c = zeros([m,size(ba,1)]);
mu_k = mu_k_SH(size(c)-1,lambda);
nu = 1.3;

% Variation Problem
% ------------------------------------------------------------------------
for i = 1:1%m(1)*m(2)
        e_j = zeros(m);
        e_j(i) = 1;
        [c,~] = newton_F(mu_k,Nop_coeff,a_bar,c,h,e_j);
%          [r, Y_0,Z_0,Z_1] = Bounds_int(intval(mu_k),intval(Nop_coeff),intval(a_bar),intval(c),intval(h),intval(e_j),intval(nu));
end