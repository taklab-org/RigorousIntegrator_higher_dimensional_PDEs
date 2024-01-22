close all
clear
clc

% Parameters 
% load equilibria1_2D_DC.mat
% N  = [6,5];
% N_Fm = [6 5 5 5 4 2];
% 
% N_bar = size(a_bar);
% a = zeros(N+1);
% a(1:6,1:6) = a_bar(:,1:6); 
% 
% a = iso_mat2vec_Fm(a,N_Fm);
% Nop_coeff = [0,0,0,-1];
% q = 2;
% L = [1,1.1];
% 
% mu_k = mu_k_dcCH_sseig(N,0.4,1,L); 
% mu_k = iso_mat2vec_Fm(mu_k,N_Fm);
% mu_k_not_in_Fm = mu_k_dcCH_sseig(N+2,0.4,1,L);
% 
% % Newton approximation 
% [a,k] = newton_F_ss(a,L,mu_k,Nop_coeff,q,N_Fm);
% 
% %Finding the eigenvalues
% 
% DF_ss = DF_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
% [V,D] = eig(DF_ss);
% [DD,I] = sort(diag(D));
% D_order = zeros(numel(mid(a)),1);
% 
% % Defining P
% nu = 1;
% omega_k = omega_weight(nu,N_Fm);
% for i = 1:numel(a)
%     vi = V(:,I(end-i+1) );
%     vi = iso_vec2mat_Fm(vi,N_Fm);
%     norm_vi = norm_1_nu_q(vi,nu,q);
%     P(:,i) = V(:,I(end-i+1)   )/(norm_vi)*omega_k(i);   
%     D_order(i) = DD(end-i+1 );
% end
% 
% % Parameters (Intervals)
% nu_int = intval(nu);
% a_int = intval(a);
% L_int = intval(L);
% mu_k_int =  iso_mat2vec_Fm_int(mu_k_SH_sseig(N,intval(0.04),L_int),N_Fm);
% Nop_coeff_int = intval(Nop_coeff);
% mu_m = max(max(abs(mu_k_int)));
% P_int = intval(P);
% 
% [Y0,Z0,Z1,Z2] = Bounds_steady_state(a_int,L_int,mu_k_int,Nop_coeff_int,nu_int,q,mu_m,N_Fm);
% if Z1 >= 1
%     error('Z_1 is greater than 1')
% end
% p = [Z2 ,-(1-Z0-Z1),Y0];
% r = sort(roots(p));
% r_int = intval(r(2));

%}
load data_DC3D.mat
% Computing extended derivative
a1 = ones(size(a));
a1 = iso_vec2mat_Fm(a1,N_Fm);
a4 = conv2(conv2(a1,a1),conv2(a1,a1));
a4(a4 > 0.9 ) = 1;
a4(a4 < 0.9 ) = 0;
N_Fm_ext = sum(a4);

a1 = ones(N(1)+1,1);
a3m = conv2(a1,conv2(a1,a1));
a3m( a3m > 0 ) = 1 ;
M3 = size(a3m)-1;
N_Fm_ext = zeros(1,M3(2)+1);

for i = 0:(size(a3m,2)-1)
    N_Fm_ext(i+1) = sum(a3m(:,i+1))-1;
end
N_Fm_ext = N_Fm_ext( N_Fm_ext > 0) ;
return

% %==============================================================
N_small = N;
N_Fm_small = N_Fm;
mu_k_small = mu_k;

 mu_inf_temp_small = flip(sort(reshape(mu_k_not_in_Fm,[],1)));
 mu_inf_small = mu_inf_temp_small(sum(N_Fm_small+1)+1);
 mu_m_small = abs(mu_inf_small );

N = [32,32];
% N_Fm = [N(1)];
% for i = 1:N(1)
%     N_Fm = [N_Fm,N(1)];
% end
% N_Fm = [30 30 30 30 30 30 30 30 29 29 28 28 28 27 26 26 25 24 23 22 21 20 19 17 16 14 11 8 3];
a_temp = iso_vec2mat_Fm(a,N_Fm_small);
N_Fm = [32 32 32 32 32 32 32 32 31 31 31 30 30 29 29 28 27 27 26 25 24 23 22 21 19 18 16 14 11 8 ];


a_small = iso_vec2mat_Fm(a,N_Fm_small);

a_temp = zeros(N+1);
a_temp(1:N_small(1)+1,1:N_small(2)+1) = a_small;
a = iso_mat2vec_Fm(a_temp,N_Fm);
mu_k = mu_k_SH_sseig(N,3,L);% lambda = 3
mu_k = iso_mat2vec_Fm(mu_k,N_Fm);

a2 = ones(N+1);
a2 = iso_vec2mat_Fm(a2,N_Fm);

a1 = ones(N_small(1)+1,1);


a3m = conv2(a2,conv2(a1,a1));
a3m( a3m > 0 ) = 1 ;
M3 = size(a3m)-1;
N_Fm_ext = zeros(1,M3(2)+1);

for i = 0:(size(a3m,2)-1)
    N_Fm_ext(i+1) = sum(a3m(:,i+1))-1;
end
N_Fm_ext = N_Fm_ext( N_Fm_ext > 0) ;



a_ext = zeros(M3(1)+1,M3(2)+1);
a_ext(1:N(1)+1,1:N(2)+1) = iso_vec2mat_Fm(a,N_Fm);
a_ext = iso_mat2vec_Fm(a_ext,N_Fm_ext);
mu_k_ext = mu_k_SH_sseig([M3(1),M3(2)],3,L);% lambda = 3
mu_k_ext = iso_mat2vec_Fm(mu_k_ext ,N_Fm_ext );

mu_k_not_in_Fm = mu_k_SH_sseig(N+2,3,L);% lambda = 3

DF_ss = DF_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
[V,D] = eig(DF_ss);
[DD,I] = sort(diag(D));
D_order = zeros(numel(mid(a)),1);
P = intval(zeros(numel(a),numel(a)));

nu = 1;

omega_k = omega_weight(nu,N_Fm);

for i = 1:numel(a)
    vi = V(:,I(end-i+1) );
    vi = iso_vec2mat_Fm(vi,N_Fm);
    norm_vi = norm_1_nu_q(vi,nu,q);
    P(:,i) = V(:,I(end-i+1)   )/(norm_vi)*omega_k(i);   
    D_order(i) = DD(end-i+1 );
end
% Raddi polynomial of Steady State
nu = intval(nu);
a = intval(a);
L = intval(L);
mu_k = intval(mu_k);
Nop_coeff = intval(Nop_coeff);
q = intval(q);

mu_inf_temp = flip(sort(reshape(mu_k_not_in_Fm,[],1)));
mu_inf = mu_inf_temp(sum(N_Fm+1)+1);
mu_m = abs(mu_inf );

% 
[Y0,Z0,Z1,Z2] = Bounds_steady_state(iso_mat2vec_Fm_int(intval(a_small),N_Fm_small),L,mu_k_small,Nop_coeff,nu,q,mu_m_small,N_Fm_small);
% 
if Z1 >= 1
    error('Z_1 is greater than 1')
end
 p = [Z2 ,-(1-Z0-Z1),Y0];
 r = sort(roots(p));


eigenvalues = diag(inv(P)*DF_ss*P);
[~,I] = sort(mid(eigenvalues));
eigenvalues = flip(eigenvalues(I));
mu_1 = eigenvalues(1);


[~,a2] = convapbqcos(intval(a_small),1,intval(a_small),1);
norm_a2 = norm_1_nu_q_int(a2,nu,q);
norm_a = norm_1_nu_q_int(intval(a_small),nu,q);
r_int = intval(r(2));
p_r = 3*(2*norm_a*r_int + r_int^2);

a_ext = intval(a_ext);
mu_k_ext = intval(mu_k_ext );
[DF11,N12,N21,~] = DF_ss_int_ext_para(a_ext,L,mu_k_ext,Nop_coeff,q,N_Fm_ext,N_Fm);

[C,lambda] = Semigroup_estimate_SH2D(a,mu_inf,nu,q,mu_k,Nop_coeff,P,N_Fm,N_Fm_ext,p_r,DF11,N12,N21);

lambda = intval(sup(lambda));
C = intval(sup(C));

r = intval(r(2));
norm_a = norm_1_nu_q_int(a,nu,q);
epsilon = intval(10^(-16));
delta = -lambda  - epsilon;

rho = inf((-3*(norm_a+r) + sqrt((-3*(norm_a+r))^2 + 4 *delta/C ))/2);

return



% N  = [10,10];
% N_Fm = [10 10 10 10 10 10 10 10 10 10 10];

N = [10,10];
N_Fm = [];
for i = 1:N(1)+1
    N_Fm = [N_Fm,N(1)];
end

N_bar = size(a_bar);

if N_bar(1) >= N(1)+1
    a = a_bar(1:N(1)+1,1:N(2)+1);
else
    a = zeros(N(1)+1,N(2)+1);
    a(1:N_bar(1),1:N_bar(2)) = a_bar;
end
a = iso_mat2vec_Fm(a,N_Fm);
Nop_coeff = [0,0,0,-1];
q = 2;
L = [1,1.1];

mu_k = mu_k_dcCH_sseig(N,0.4,1,L); % lambda = 3
mu_k = iso_mat2vec_Fm(mu_k,N_Fm);
mu_k_not_in_Fm = mu_k_dcCH_sseig(N+2,0.4,1,L);% lambda = 3


[a,k] = newton_F_ss(a,L,mu_k,Nop_coeff,q,N_Fm);

% %==============================================================
N_small = N;
N_Fm_small = N_Fm;
mu_k_small = mu_k;

 mu_inf_temp_small = flip(sort(reshape(mu_k_not_in_Fm,[],1)));
 mu_inf_small = mu_inf_temp_small(sum(N_Fm_small+1)+1);
 mu_m_small = abs(mu_inf_small );


% N_Fm = [N(1)];
% for i = 1:N(1)
%     N_Fm = [N_Fm,N(1)];
% end
% N_Fm = [30 30 30 30 30 30 30 30 29 29 28 28 28 27 26 26 25 24 23 22 21 20 19 17 16 14 11 8 3];
a_temp = iso_vec2mat_Fm(a,N_Fm_small);



a_small = iso_vec2mat_Fm(a,N_Fm_small);

a_temp = zeros(N+1);
a_temp(1:N_small(1)+1,1:N_small(2)+1) = a_small;
a = iso_mat2vec_Fm(a_temp,N_Fm);
mu_k = mu_k_dcCH_sseig(N,0.4,1,L);% lambda = 3
mu_k = iso_mat2vec_Fm(mu_k,N_Fm);

a2 = ones(N+1);
a2 = iso_vec2mat_Fm(a2,N_Fm);

a1 = ones(N_small(1)+1,N_small(1)+1);


a3m = conv2(a2,conv2(a1,a1));
a3m( a3m > 0 ) = 1 ;
M3 = size(a3m)-1;

N_ext = N + [20,20];
N_Fm_ext = [];
for i = 1:N_ext(1)+1
    N_Fm_ext = [N_Fm_ext,N_ext(1)];
end



% N_Fm_ext = zeros(1,M3(2)+1);
% for i = 0:(size(a3m,2)-1)
%     N_Fm_ext(i+1) = sum(a3m(:,i+1))-1;
% end
% N_Fm_ext = N_Fm_ext( N_Fm_ext > 0) ;



a_ext = zeros(M3(1)+1,M3(2)+1);
a_ext(1:N(1)+1,1:N(2)+1) = iso_vec2mat_Fm(a,N_Fm);
a_ext = iso_mat2vec_Fm(a_ext,N_Fm_ext);
mu_k_ext = mu_k_SH_sseig([M3(1),M3(2)],3,L);% lambda = 3
mu_k_ext = iso_mat2vec_Fm(mu_k_ext ,N_Fm_ext );





mu_k_not_in_Fm = mu_k_SH_sseig(N+2,3,L);% lambda = 3

DF_ss = DF_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
[V,D] = eig(DF_ss);
[DD,I] = sort(diag(D));
D_order = zeros(numel(mid(a)),1);
P = intval(zeros(numel(a),numel(a)));

nu = 1;

omega_k = omega_weight(nu,N_Fm,q);

for i = 1:numel(a)
    vi = V(:,I(end-i+1) );
    vi = iso_vec2mat_Fm(vi,N_Fm);
    norm_vi = norm_1_nu_q(vi,nu,q);
    P(:,i) = V(:,I(end-i+1)   );%/(norm_vi)*omega_k(i);   
    D_order(i) = DD(end-i+1 );
end
% Raddi polynomial of Steady State
nu = intval(nu);
a = intval(a);
L = intval(L);
mu_k = intval(mu_k);
Nop_coeff = intval(Nop_coeff);
q = intval(q);

mu_inf_temp = flip(sort(reshape(mu_k_not_in_Fm,[],1)));
mu_inf = mu_inf_temp(sum(N_Fm+1)+1);
mu_m = abs(mu_inf );

% % 
% [Y0,Z0,Z1,Z2] = Bounds_steady_state(iso_mat2vec_Fm_int(intval(a_small),N_Fm_small),L,mu_k_small,Nop_coeff,nu,q,mu_m_small,N_Fm_small);
% % 
% if Z1 >= 1
%     error('Z_1 is greater than 1')
% end
%  p = [Z2 ,-(1-Z0-Z1),Y0];
%  r = sort(roots(p));
r = 10^-5;

% Finding C1 and C2

DF_ss_int = DF_steady_state_int(a,L,mu_k,Nop_coeff,q,N_Fm);
P_int = intval(P);
P_inv = inv(P_int);
C1 = max([norm_B_nu_q_int(P,nu,N,q,N_Fm),1]);
C2 = max([norm_B_nu_q_int(P_inv,nu,N,0,N_Fm),1]);


% Finding <n>^q and lambda^_n for fixed t
eigenvalues = diag(P_inv*DF_ss_int*P_int);
lambda_hat_k = intval(abs(iso_vec2mat_Fm(mu_k_ext,N_Fm_ext)));
lambda_hat_k(1:N(1)+1,1:N(2)+1) = abs(iso_vec2mat_Fm_int(eigenvalues,N_Fm));
% [brac_n,lambda_hat_n] = constant_eta(iso_vec2mat_Fm_int(lambda_hat_k,N_Fm_ext),t,q);


% Computing || E || 
[~,a2] = convapbqcos(intval(a_small),1,intval(a_small),1);
norm_a2 = norm_1_nu_q_int(a2,nu,0);
norm_a = norm_1_nu_q_int(intval(a_small),nu,0);
r_int = intval(r); %intval(r(2));
p_r = 3*(2*norm_a*r_int + r_int^2);


norm_E = E_cally(a_ext,L,mu_k_ext,Nop_coeff,q,N_Fm_ext,N_Fm,nu,p_r,P_int,P_inv,a,mu_k);


return
[~,I] = sort(mid(eigenvalues));
eigenvalues = flip(eigenvalues(I));
mu_1 = eigenvalues(1);




a_ext = intval(a_ext);
mu_k_ext = intval(mu_k_ext );
[DF11,N12,N21,~] = DF_ss_int_ext_para(a_ext,L,mu_k_ext,Nop_coeff,q,N_Fm_ext,N_Fm);

[C,lambda] = Semigroup_estimate_DC2D(a,mu_inf,nu,q,mu_k,Nop_coeff,P,N_Fm,N_Fm_ext,p_r,DF11,N12,N21);

% lambda = intval(sup(lambda));
% C = intval(sup(C));
% 
% r = intval(r(2));
% norm_a = norm_1_nu_q_int(a,nu,q);
% epsilon = intval(10^(-16));
% delta = -lambda  - epsilon;
% 
% rho = inf((-3*(norm_a+r) + sqrt((-3*(norm_a+r))^2 + 4 *delta/C ))/2);
% a = iso_vec2mat_Fm_int(a_small,N_Fm_small);




