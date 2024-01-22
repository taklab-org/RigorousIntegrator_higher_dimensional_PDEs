close all
clear
clc


% load equilibria1_2D_SH.mat
% load equilibria2_2D_SH.mat
load equilibria_3D_SH.mat
q = 0;
% N  = [9,9];
N  = [3,3,3];
a = zeros(N(1)+1,N(2)+1,N(3)+1);
N_bar = size(a_bar);
a(1:N_bar(1),1:N_bar(2),1:N_bar(3)) = a_bar;
L = [1,1.1,1.2];

mu_k = mu_k_SH_sseig(N,0.04,L);

mu_m = abs(mu_k(end,1));
Nop_coeff = [0,0,0,-1];

fprintf("\nFinding Equilibrium\n")
[a,~] = newton_F_ss(a,L,mu_k,Nop_coeff,q);

% plot_SH_profile(mid(a),[1,1.1]),pause(0.01)

% Computation of the eigenfunctions
DF_ss = DF_steady_state(a,L,mu_k,Nop_coeff,q);
[V,D] = eig(DF_ss); 
[DD,I] = sort(diag(D));
D_order = zeros(size(a));
P = intval(zeros(numel(a),numel(a)));
for i = 1:numel(a)
    P(:,i) = V(:,I(end-i+1))/(sum(abs(V(:,I(end-i+1)))));   
    D_order(i) = D(I(end-i+1),I(end-i+1));
end

% Raddi polynomial of Steady State
nu = intval(1.0);
a = intval(a);
L = intval(L);
mu_k = intval(mu_k);
Nop_coeff = intval(Nop_coeff);
q = intval(q);
mu_m = intval(mu_m);

tic
[Y0,Z0,Z1,Z2] = Bounds_steady_state(a,L,mu_k,Nop_coeff,nu,q,mu_m);
toc

if Z1 >= 1
    error('Z_1 is greater than 1')
end
p = [Z2 ,-(1-Z0-Z1),Y0];
r = sort(roots(p));



tic
% [G_0,C,rho,epsilon,delta]  =Gershgorin_General(a,P,mu_k,Nop_coeff,L,intval(r(2)),nu,mu_k(end,1),q);
 [G_0,C,rho,epsilon,delta]  = Gershgorin_ring_3D_SH(a,P,mu_k,Nop_coeff,L,intval(r(2)),nu,mu_k(end,1),q);
toc





