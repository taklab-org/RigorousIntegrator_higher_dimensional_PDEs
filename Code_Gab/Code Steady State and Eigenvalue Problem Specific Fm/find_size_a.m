close all
clear
clc

load equilibria2_2D_DC.mat
% load equilibria2_2D_SH.mat
%load equilibria2_2D_SH_alt.mat
load Map_Fm.mat
Map_Fm = mid(Map_Fm);

N  = [5,5];
N_Fm = [5 5 5 5 5 5];

% N  = [18,18]; 
% 
% N_Fm = [18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18];
% a = zeros(N(1)+1,N(2)+1);
N_bar = size(a_bar);

if N_bar(1) >= N(1)+1
    a = a_bar(1:N(1)+1,1:N(2)+1);
else
    a = zeros(N(1)+1,N(2)+1);
    a(1:N_bar(1),1:N_bar(2)) = a_bar;
end
a = iso_mat2vec_Fm(a,N_Fm);


L  = [1,1.1];
mu_k = mu_k_dcCH_sseig(N,0.4,1,L);
mu_k = iso_mat2vec_Fm(mu_k,N_Fm);

mu_k_not_in_Fm = mu_k_dcCH_sseig(N+2,0.4,1,L);


Nop_coeff = [0,0,0,-1];
q = 2;

[F] = F_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
DF =  DF_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
DF_fin = fin_dif_ss(a,L,mu_k,Nop_coeff,q,N_Fm);


[a,k] = newton_F_ss(a,L,mu_k,Nop_coeff,q,N_Fm);


N_small = N;
N_Fm_small = N_Fm;

N = [6,6];
N_Fm = [N(1)];
for i = 1:N(1)
    N_Fm = [N_Fm,N(1)];
end
%  N_Fm = [30 30 30 30 30 30 30 30 29 29 28 28 28 27 26 26 25 24 23 22 21 20 19 17 16 14 11 8 3];

a_temp = iso_vec2mat_Fm(a,N_Fm_small);

a = zeros(N+1);

a(1:size(a_temp,1) , 1:size(a_temp,1) ) = a_temp;

a = iso_mat2vec_Fm(a,N_Fm);
mu_k = mu_k_dcCH_sseig(N,0.4,1,L);
mu_k = iso_mat2vec_Fm(mu_k,N_Fm);
mu_k_not_in_Fm = mu_k_SH_sseig(N+2,3,L);

DF_ss = DF_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
[V,D] = eig(DF_ss);
[DD,I] = sort(diag(D));
D_order = zeros(numel(mid(a)),1);
P = intval(zeros(numel(a),numel(a)));

for i = 1:numel(a)
    vi = V(:,I(end-i+1) );
    vi = iso_vec2mat_Fm(vi,N_Fm);
    norm_vi = norm_1_nu_q(vi,1,q);
    P(:,i) = V(:,I(end-i+1)   )/norm_vi;   
    D_order(i) = DD(end-i+1 );
end
diff = 0;
sort(abs(mu_k));
for i = 1:length(D_order)-1
    delta = abs(abs(D_order(i)) - D_order(i+1));
    if delta > diff
    diff = delta;
    k_ind = i;
    end
end
mu_inf_temp = flip(sort(reshape(mu_k_not_in_Fm,[],1)));
mu_inf = mu_inf_temp(sum(N_Fm+1)+1);
mu_m = abs(mu_inf );

A = reshape(mu_k_not_in_Fm,N(1)+3,N(2)+3);
mu_inf = A(N(1)+2,1);
 
eigenvalue = D_order( D_order > mu_inf) ;

epsilon = sum(1./(abs(mu_inf)-40 - abs(eigenvalue)))

1/epsilon

D_size = reshape(mu_k,N+1);
D_size(D_size < mu_inf ) = 0;
D_size(D_size ~= 0  ) = 1;

N_Fm = zeros(1,size(D_size,2));
for i = 1:size(D_size,2)
    N_Fm(i) = sum(D_size(:,i))-1;
end
N_Fm = N_Fm(N_Fm ~= -1);