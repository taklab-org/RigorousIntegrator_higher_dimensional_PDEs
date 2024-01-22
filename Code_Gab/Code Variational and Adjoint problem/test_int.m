clear
close all
clc

%==========================================================================
% 2D SH
%==========================================================================
load('data_gab_rand')
% a_bar = a_bar(1:3,1:3,:);
N = size(a_bar)-1;
Nop_coeff = [0,0,0,-1];
% h = 0.1;
unit_pad = zeros(size(N));
unit_pad(end) = 7;
a_bar = padarray(a_bar,unit_pad,'post');
c = a_bar(1:3,1:3,:);

N = size(a_bar)-1;
Phi = zeros(numel(c),prod(N(1:end-1)));
mu_k = mu_k_SH(N,lambda);

% r = zeros(1,prod(N(1:end-1)));
% Y_0 = zeros(1,prod(N(1:end-1)));
% Z_0 = zeros(1,prod(N(1:end-1)));
% Z_1 = zeros(1,prod(N(1:end-1)));
nu = 1.5;


for i = 1:1%((N(1)+1)*(N(2)+1))
        e_j = zeros(N(1)+1,N(2)+1);
        e_j(i) = 1;
        [c,~] = newton_F(mu_k,Nop_coeff,a_bar,c,h,e_j);
        [r, Y0,Z0,Z1] = rigorous_bounds_lin_SH(c,a_bar,h,e_j,lambda,nu);
        return
        Phi(:,i) = reshape(c,[numel(c),1]);
        F1 = F_lin_fin_dim(mu_k,Nop_coeff,a_bar,c,h,e_j);
        F2 = F_lin_fin_dim_int(intval(mu_k),intval(Nop_coeff),intval(a_bar),intval(c),intval(h),intval(e_j));
        F3 = F_lin_fin_dim_int2(intval(mu_k),intval(Nop_coeff),intval(a_bar),intval(c),intval(h),intval(e_j));
        [r(i), Y_0(i),Z_0(i),Z_1(i)] = Bounds_int(intval(mu_k),intval(Nop_coeff),intval(a_bar),intval(c),intval(h),intval(e_j),intval(nu));
end
return
%==========================================================================
% 3D SH
%==========================================================================
% load('gab_data_test_3D.mat')
% Nop_coeff = [0,0,0,-1];
% a_bar = permute(ba,[2,3,4,1]);
% % a_bar = a_bar(1:3,1:3,1:3,1:10);
% N = size(a_bar)-1;
% mu_k = mu_k_SH(N,lambda);  
% e_j = zeros(N(1)+1,N(2)+1,N(3)+1);
% e_j(1)= 1;
% nu = 1.5;
% [c,~] = newton_F(mu_k,Nop_coeff,a_bar,a_bar,h,e_j);
% [r, Y0,Z0,Z1] = rigorous_bounds_lin_SH(c,a_bar,h,e_j,lambda,nu);





