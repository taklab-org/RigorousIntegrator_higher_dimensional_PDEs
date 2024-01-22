clc
close all
clear


%test 2d Sh
load equilibria2_2D_SH.mat
N  = [16,16];
a = zeros(N(1)+1,N(2)+1);
N_bar = size(a_bar);
a(1:N_bar(1),1:N_bar(2)) = a_bar;
L = [1,1];
mu_k = mu_k_SH_sseig(N,3,L);
Nop_coeff = [0,0,0,-1];

F_ss = F_steady_state(a,L,mu_k,Nop_coeff);




[a,~] = newton_F_ss(a,L,mu_k,Nop_coeff);
norm(a)

DF_ss = reshape(DF_steady_state(a,L,mu_k,Nop_coeff),numel(a),numel(a));
DF_fin_ss = reshape(fin_dif_ss(a,L,mu_k,Nop_coeff),numel(a),numel(a));

[V,D] = eig(DF_ss);
[lambda,k] = max(diag(D));
phi = reshape(V(:,k),N+1); phi(1) = 0;

[F] = F_eigenpair(a,phi,lambda,L,mu_k,Nop_coeff);
[phi,lambda,k] = newton_F_eig(a,phi,lambda,L,mu_k,Nop_coeff);

F_all = F_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff);

x = [a,phi];
DF_all = reshape(DF_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff),numel(x),numel(x));
[a,phi,lambda,k] = newton_F_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff);

% phi = rand(size(a)); 
% lambda = -1;
% 
% [phi,lambda,k] = newton_F_eig(a,phi,lambda,L,mu_k,Nop_coeff);



% %test 3d Sh
% load equilibria_3D_SH.mat
% N  = [3,3,3];
% a = rand(N+1);
% L = [1,1,1];
% mu_k = mu_k_SH_sseig(N,3,L);
% 
% DF = reshape(DF_steady_state(a,L,mu_k,Nop_coeff),numel(a),numel(a));
% DF_fin = reshape(fin_dif_ss(a,L,mu_k,Nop_coeff),numel(a),numel(a));
% norm(DF -DF_fin )
% surf(abs(DF -DF_fin))
% 
% return
% N_bar = size(a_bar);
% a(1:N_bar(1),1:N_bar(2),1:N_bar(3)) = a_bar;
% L = [1,1.1,1.2];
% mu_k = mu_k_SH_sseig(N,3,L);
% Nop_coeff = [0,0,0,-1];
% nop = N_op_nonlin(Nop_coeff,a);
% 
% F_ss = F_steady_state(a,L,mu_k,Nop_coeff);
% [a,~] = newton_F_ss(a,L,mu_k,Nop_coeff);
% norm(reshape(a,numel(a),1))
% 
% DF = reshape(DF_steady_state(a,L,mu_k,Nop_coeff),numel(a),numel(a));
% DF_fin = reshape(fin_dif_ss(a,L,mu_k,Nop_coeff),numel(a),numel(a));
% 







