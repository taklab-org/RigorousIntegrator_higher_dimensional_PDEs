clc
close all
clear

load equilibria_3D_SH_more_modes.mat
N  = [6,6,6];
a = zeros(N(1)+1,N(2)+1,N(3)+1);
N_bar = size(a_bar);
a(1:N_bar(1),1:N_bar(2),1:N_bar(3)) = a_bar;
L = [1,1.1,1.2];
mu_k = mu_k_SH_sseig(N,0.4,L);
Nop_coeff = [0,0,0,-1];

[u,D_2u,D_4u]  = plot_cos_3d(a,L);
V = @(lambda) (lambda-1)*u(1,1,1) - 2*D_2u(1,1,1) - D_4u(1,1,1) - u(1,1,1).^3;
