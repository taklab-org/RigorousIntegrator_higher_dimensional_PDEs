
close all
clear
clc


%load equilibria1_2D_SH.mat
load equilibria2_2D_SH.mat
%load equilibria2_2D_SH_alt.mat
q = 0;
% N  = [9,9];
N  = [9,9];
% a = zeros(N(1)+1,N(2)+1);
N_bar = size(a_bar);


if N_bar(1) >= N(1)+1
    a = a_bar(1:N(1)+1,1:N(2)+1);
else
    a = zeros(N(1)+1,N(2)+1);
    a(1:N_bar(1),1:N_bar(2)) = a_bar;

end


L = [1,1.1];
mu_k = mu_k_SH_sseig(N,5,L);

mu_k_not_in_Fm = mu_k_SH_sseig(N+1,5,L);



mu_m = abs(mu_k(end,1));
Nop_coeff = [0,0,0,-1];

fprintf("\nFinding Equilibrium\n")
[a,~] = newton_F_ss(a,L,mu_k,Nop_coeff,q,N_Fm);

% plot_SH_profile(mid(a),[1,1.1]),pause(0.01)

% Computation of the eigenfunctions
DF_ss = DF_steady_state(a,L,mu_k,Nop_coeff,q);
[V,D] = eig(DF_ss);
[DD,I] = sort(diag(D));
D_order = zeros(numel(mid(a)),1);
P = intval(zeros(numel(a),numel(a)));


M_ind = zeros(N+1);
k_id = 0;
Vect_list_idx = (1:numel(a))';
N_delta = abs(N(1) - N(2)) +1;
for i = 1:min(N +1)

    if i ~= min(N +1)
        v1 = Vect_list_idx( 1: 1 + k_id);
        v2 = Vect_list_idx( end - k_id: end);
        k_id = k_id + 1;
        Vect_list_idx = Vect_list_idx(1+k_id:end-k_id);
        M_ind = M_ind + flip(diag(v1, -min(N +1 ) + i)) + flip(diag(v2, min(N +1 ) - i));
    else
        for j = 1:N_delta
            v1 = Vect_list_idx( 1: 1 + k_id);
            Vect_list_idx = Vect_list_idx(1+k_id:end);
            M_ind = M_ind + flip(diag(v1, -min(N +1 ) + i)) ;
        end
    end
end

M_ind  = reshape(M_ind,[],1);

DD = flip(DD);
for i = 1:numel(a)
    vi = V(:,M_ind(i) );
    vi = reshape(vi,N+1);
    norm_vi = norm_1_nu_q(vi,1,q);
    P(:,i) = V(:,M_ind(i)  )/norm_vi;   
    D_order(i) = DD(M_ind(i));

    % vi = V(:,I(end-i+1));
    % vi = reshape(vi,N+1);
    % norm_vi = norm_1_nu_q(vi,1,q);
    % P(:,i) = V(:,I(end-i+1))/norm_vi;   
    % D_order(i) = D(I(end-i+1),I(end-i+1));
end
% Raddi polynomial of Steady State
nu = intval(1);
a = intval(a);
L = intval(L);
mu_k = intval(mu_k);
Nop_coeff = intval(Nop_coeff);
q = intval(q);
mu_m = intval(mu_m);
% tic
% [Y0,Z0,Z1,Z2] = Bounds_steady_state(a,L,mu_k,Nop_coeff,intval(1),q,mu_m);
% toc
% if Z1 >= 1
%     error('Z_1 is greater than 1')
% end
% p = [Z2 ,-(1-Z0-Z1),Y0];
% r = sort(roots(p));

r = [1,0.25813572708351*10^(-8)];
%
% tic
%  [G_0,C,rho,epsilon,delta]  = Gershgorin_ring(a,P,mu_k,Nop_coeff,L,intval(r(2)),nu,mu_k(end,1),q);
%  G_0_test_1 = G_0;
%  toc
%  tic
% [G_0,C,rho,epsilon,delta]  =Gershgorin_General(a,P,mu_k,Nop_coeff,L,intval(r(2)),nu,mu_k(end,1),q);
% G_0_test_2 = G_0;
% toc

mu_inf_temp = flip(sort(reshape(mu_k_not_in_Fm,[],1)));
mu_inf = mu_inf_temp(81);



eigenvalues = diag(inv(P)*DF_ss*P);
[~,I] = sort(mid(eigenvalues));
eigenvalues = flip(eigenvalues(I));
eigenvalues = eigenvalues(1:80);
mu_1 = eigenvalues(1);

[C] = Semigroup_estimate_SH2D(a,mu_1,mu_inf,nu,q,L,mu_k,Nop_coeff,P,eigenvalues)


