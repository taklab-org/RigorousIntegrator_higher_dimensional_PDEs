close all
clear
clc

%Loading Data
tic
load equilibria_3D_SH.mat

%Parameters
N  = [3,2,2];
L = [1,1.1,1.2];
mu_k = mu_k_SH_sseig(N,0.04,L);
mu_index = abs(mu_k(end,1,1));
N_Fm_map = mu_k;
N_Fm_map(abs(N_Fm_map) <= mu_index) = 1;
N_Fm_map(abs(N_Fm_map) > mu_index) = 0;


N_bar = size(a_bar);
a = zeros(size(N_Fm_map ));
a(1:N_bar(1),1:N_bar(2),1:N_bar(3)) = a_bar;
a = iso_3D_mat2vec_Fm(a,N_Fm_map);
Nop_coeff = [0,0,0,-1];
q = 0;
mu_k = mu_k_SH_sseig(N,0.04,L); % lambda = 3
mu_k = iso_3D_mat2vec_Fm(mu_k,N_Fm_map);
mu_k_not_in_Fm = mu_k_SH_sseig(N+2,0.04,L);% lambda = 3

% Newton to converge to better approximation of solution
[a,k] = newton_F_ss(a,L,mu_k,Nop_coeff,q,N_Fm_map);

%Finding eigenvalues
DF_ss = DF_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm_map);
[V,D] = eig(DF_ss);
[DD,I] = sort(diag(D));
D_order = zeros(numel(mid(a)),1);

% Defining P
nu = 1;
q = 0;
omega_k = omega_weight_3D(nu,N_Fm_map,q);
for i = 1:numel(a)
    vi = V(:,I(end-i+1) );
    norm_vi = norm_1_nu_q_3D(vi,nu,N_Fm_map,q);
    P(:,i) = V(:,I(end-i+1)   )/(norm_vi)*omega_k(i);   
    D_order(i) = DD(end-i+1 );
end

% Parameters (Intervals)
nu_int = intval(nu);
a_int = intval(a);
L_int = intval(L);
mu_k_int =  iso_3D_mat2vec_Fm(mu_k_SH_sseig(N,intval(0.04),L_int),N_Fm_map);
Nop_coeff_int = intval(Nop_coeff);
mu_m = max(max(max(abs(mu_k_int))));
P_int = intval(P);


% Proving steady state of finite truncation in time
[Y0,Z0,Z1,Z2] = Bounds_steady_state_3D(a_int,L_int,mu_k_int,Nop_coeff_int,nu_int,q,mu_m,N_Fm_map);
if Z1 >= 1
    error('Z_1 is greater than 1')
end
p = [Z2 ,-(1-Z0-Z1),Y0];
r = sort(roots(p));
r_int = intval(r(2));

% Computing extended derivative
one_a = N_Fm_map; 
[~,one_a4] = convapbqcos(one_a,2,one_a,2);
one_a4(one_a4 >= 0.9) = 1;
one_a4(one_a4 < 0.9) = 0;
N_Fm_map_ext = one_a4;

a_mat = iso_3D_vec2mat_Fm(a_int,N_Fm_map);
a_pad = intval(zeros(size(N_Fm_map_ext)));
a_pad(1:N(1)+1,1:N(2)+1,1:N(3)+1) = a_mat;
a_pad = iso_3D_mat2vec_Fm(a_pad,N_Fm_map_ext);

mu_k_ext = mu_k_SH_sseig(size(N_Fm_map_ext)-1,0.04,L_int);
mu_k_ext = iso_3D_mat2vec_Fm(mu_k_ext,N_Fm_map_ext);

[D11,D12,D21] = DF_split_3D_V2(a_pad,L_int,mu_k_ext,Nop_coeff,q,N_Fm_map_ext,N_Fm_map);

% Finding the residu polynomial
mu_k_ext = iso_3D_vec2mat_Fm(mu_k_ext,N_Fm_map_ext);
N_Fm_pad = zeros(size(N_Fm_map_ext));
Size_N_Fm = size(N_Fm_map);
N_Fm_pad(1:Size_N_Fm(1),1:Size_N_Fm(2),1:Size_N_Fm(3)) = N_Fm_map;
mu_k_delta = (ones(size(N_Fm_pad)) - N_Fm_pad ).*mu_k_ext;
mu_k_delta(mu_k_delta >= 0) = NaN;
mu_inf = min(min(min(abs(mu_k_delta))));
norm_a_bar = norm_1_nu_q_3D(a_int,nu_int,N_Fm_map,0);
Norm_DN = 3*norm_a_bar^2+4*r_int*norm_a_bar+2*r_int^2;
p_r = r_int*norm_a_bar+2*r_int^2;

% Semigroup estimate
[C,lambda] = Semigroup_estimate_SH3D(D11,D12,D21,P_int,mu_inf,Norm_DN,nu_int,N_Fm_map,N_Fm_map_ext,p_r);

%Computing  rho 
epsilon = intval(10^(-16));
delta = lambda  - epsilon;
rho = inf((-3*(norm_a_bar+r_int) + sqrt((3*(norm_a_bar+r_int))^2 + 4 *delta/C ))/2);
toc