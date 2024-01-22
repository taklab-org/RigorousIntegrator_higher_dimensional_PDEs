clear
close all
clc

profile off
profile on

%=================================================================
% 3D SH
%=================================================================
% load test_SH_3d
% 
% a = inf(a);
% 
% ell = 11;
% m1 = 5;
% m2 = 5;
% m3 = 5;
% M = [m1,m2,m3,ell];
% [n,m,k,ELL] = ndgrid(0:m1,0:m2,0:m3,0:ell);
% c = rand(m1+1,m2+1,m3+1,ell+1)./(4.^(n+m+k+ELL));
% 
% N1 = 11;
% N2 = 11;
% N3 = 11;
% N = [N1,N2,N3,ell];
% [n,m,k,ELL] = ndgrid(0:N1,0:N2,0:N3,0:ell);
% a = rand(N1+1,N2+1,N3+1,ell+1)./(4.^(n+m+k+ELL));
% 
% Nop_coeff = [0,0,0,-1];
% lambda = 1;
% mu_k= mu_k_SH(M,lambda);
% a_bar = a;
% h = 2^-8;
% e_j = zeros(N1+1,N2+1,N3+1); e_j(1) = 1;
% 
% clearvars -except mu_k Nop_coeff a_bar c h e_j
% [F] = F_lin_fin_dim(mu_k,Nop_coeff,a_bar,c,h,e_j);
% 
% profile viewer
%=================================================================
% 2D SH
%=================================================================

% load 'test_2d.mat'
load ba_2D_SH.mat
a = permute(ba,[2,3,1]);
N = N';
c = zeros(m');
lambda = 3;
% a = a(1:N1+1,1:N2+1,1:N3+1);
% N = [N1,N2,N3];
% [n,m,ELL] = ndgrid(0:N(1),0:N(2),0:N(3));
% c = rand(N(1)+1,N(2)+1,N(3)+1)./(4.^(n+m+ELL));
Nop_coeff = [0,0,0,-1];
mu_k = mu_k_SH(N,lambda);

a_bar = a;
h =  2^-5;

e_j = zeros(N(1)+1,N(2)+1); e_j(1) = 1;
F = F_lin_fin_dim(mu_k,Nop_coeff,a_bar,c,h,e_j);
[c,k] = newton_F(mu_k,Nop_coeff,a_bar,c,h,e_j);
save('DATA_2D_SH','mu_k','Nop_coeff','a_bar','c','h','e_j')


nu = 1.1;
[r, Y_0,Z_0,Z_1] = Bounds(mu_k,Nop_coeff,a_bar,c,h,e_j,nu);



