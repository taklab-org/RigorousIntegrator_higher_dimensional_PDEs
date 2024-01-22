close all
clear
clc

%Loading Data
load equilibria1_2D_SH.mat
tic
%Parameters
N  = [9,9];
L = [1,1.1];
mu_k = mu_k_SH_sseig(N,3,L);
mu_index = abs(mu_k(end,1));
N_Fm_map = [9 9 9 9 8 8 7 6 4 1];

N_bar = size(a_bar);
a = zeros(N_Fm_map(1)+1,size(N_Fm_map,2));
a(1:N_bar(1),1:N_bar(2)) = a_bar;
a = iso_mat2vec_Fm(a,N_Fm_map);
Nop_coeff = [0,0,0,-1];
q = 0;
mu_k = mu_k_SH_sseig(N,3,L); % lambda = 3
mu_k = iso_mat2vec_Fm(mu_k,N_Fm_map);
mu_k_not_in_Fm = mu_k_SH_sseig(N+2,3,L);% lambda = 3

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
omega_k = omega_weight(nu,N_Fm_map,q);
for i = 1:numel(a)
    vi = iso_vec2mat_Fm(V(:,I(end-i+1) ), N_Fm_map );
    norm_vi = norm_1_nu_q(vi,nu,q);
    P(:,i) = V(:,I(end-i+1)   )/(norm_vi)*omega_k(i);   
    D_order(i) = DD(end-i+1 );
end

% Parameters (Intervals)
nu_int = intval(nu);
a_int = intval(a);
L_int = intval(L);
mu_k_int =  iso_mat2vec_Fm_int(mu_k_SH_sseig(N,intval(3),L_int),N_Fm_map);
Nop_coeff_int = intval(Nop_coeff);
mu_m = max(max(max(abs(mu_k_int))));
P_int = intval(P);

% Proving steady state of finite truncation in time
[Y0,Z0,Z1,Z2] = Bounds_steady_state(a_int,L_int,mu_k_int,Nop_coeff_int,nu_int,q,mu_m,N_Fm_map);
if Z1 >= 1
    error('Z_1 is greater than 1')
end
p = [Z2 ,-(1-Z0-Z1),Y0];
r = sort(roots(p));
r_int = intval(r(2));


% Computing extended derivative
ones_mat = zeros(N+1);
for i = 1:length(N_Fm_map)
    for j = 0:N_Fm_map(i)
        ones_mat(j+1,i) = 1 ;
    end
end

[~,ones_4] = convapbqcos(ones_mat,2,ones_mat,2); 
ones_4(ones_4 >= 0.9) = 1;
ones_4(ones_4 < 0.9) = 0;

N_Fm_map_ext = sum(ones_4)-1;
a_pad = intval(zeros(size(ones_4)));
a_mat = iso_vec2mat_Fm_int(a_int,N_Fm_map);
a_pad(1:N(1)+1,1:N(2)+1) = a_mat;
a_pad = iso_mat2vec_Fm_int(a_pad,N_Fm_map_ext);
N_ext = size(ones_4)-1;

mu_k_ext = mu_k_SH_sseig(N_ext,3,L_int); % lambda = 3
mu_k_ext = iso_mat2vec_Fm_int(mu_k_ext,N_Fm_map_ext);

[D11,D12,D21] = DF_split_2D(a_pad,L_int,mu_k_ext,Nop_coeff,q,N_Fm_map_ext,N_Fm_map);

% Semigroup estimate
norm_a_bar = norm_1_nu_q_int(a_mat,nu_int,0);
p_r = r_int*norm_a_bar+2*r_int^2;
Norm_DN = 3*norm_a_bar^2+4*r_int*norm_a_bar+2*r_int^2;
%--------------------------------------------------------------------------
N_Fm_maping = zeros(N_Fm_map(1)+1,length(N_Fm_map));
N_Fm_ext_map = zeros(N_Fm_map_ext(1)+1,length(N_Fm_map_ext));

for i = 1:length(N_Fm_map)
    for j = 0:N_Fm_map(i)
        N_Fm_maping(j+1,i) = 1 ;
    end
end
for i = 1:length(N_Fm_map_ext)
    for j = 0:N_Fm_map_ext(i)
        N_Fm_ext_map(j+1,i) = 1 ;
    end
end

N_Fm_map_pad = zeros(size(N_Fm_ext_map));
N_Fm_size = size(N_Fm_maping);
N_Fm_map_pad(1:N_Fm_size(1),1:N_Fm_size(2)) = N_Fm_maping;

N_Fm_delta = iso_mat2vec_Fm( N_Fm_ext_map - N_Fm_map_pad,N_Fm_map_ext) ;
mu_k_delta = mu_k_ext.*N_Fm_delta ;
mu_k_delta(mu_k_delta >= 0) = NaN;
mu_inf = min(abs(mu_k_delta));
%--------------------------------------------------------------------------
[C,lambda] = Semigroup_estimate_SH2D_v2(D11,D12,D21,P_int,mu_inf,Norm_DN,nu_int,N_Fm_map,N_Fm_map_ext,p_r);

%Computing  rho 
epsilon = intval(10^(-16));
delta = lambda  - epsilon;
rho = inf((-3*(norm_a_bar+r_int) + sqrt((3*(norm_a_bar+r_int))^2 + 4 *delta/C ))/2);
toc