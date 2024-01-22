function [Y0,Z0,Z1,Z2] = Bounds_steady_state_3D(a,L,mu_k,Nop_coeff,nu,q,mu_m,N_Fm)
a_vect = a;
mu_k_vec = mu_k;

a = iso_3D_vec2mat_Fm(a,N_Fm);
mu_k = iso_3D_vec2mat_Fm(mu_k,N_Fm);

% Parameters
Na_num = numel(mid(a));
Na = size(a)-1;
Ma = size(a);
for i = 1:length(Ma)
    otherdims_a{i} = 1:Ma(i);
end

Ma = (size(Nop_coeff,2)-1)*Na;
a_pad = intval(zeros(Ma+1));
a_pad(otherdims_a{:}) = a;


% Ma = (size(Nop_coeff,2)-1)*Na; % Size of the extension(e.g. Ma = 3*Na for 2D SH)
% unit_pad = (Ma-Na).*ones(1,length(size(a))); %use to pad a
% a_pad = padarray(a,unit_pad ,'post');% Padded a_bar 
% phi_pad = padarray(phi,unit_pad ,'post');% Padded phi_bar 

% Operators

F1N = F_steady_state_int(a_vect,L,mu_k,Nop_coeff,q,N_Fm);
DFN = DF_steady_state_int(mu_k_vec,L,mu_k,Nop_coeff,q,N_Fm);

AN = inv(mid(DFN));


nop_ext = N_op_nonlin_int(Nop_coeff,a_pad);
nop_ext = banach_alg(a_pad,Nop_coeff,q,nop_ext);
M_ext = size(nop_ext)-1;
dim = length(M_ext);
s = arrayfun(@(k) 0:k, M_ext(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if q >0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^q;
    end
else
    K_sq = intval(ones(size(k{1})));
end

nop_ext = K_sq.*nop_ext;

    






nop= N_op_nonlin_int(Nop_coeff,a);
M = size(nop)-1;
dim = length(M);
s = arrayfun(@(k) 0:k, M(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if q >0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^q;
    end
else
    K_sq = intval(ones(size(k{1})));
end
nop = K_sq.*nop;

nop_pad = intval(zeros(Ma+1));
nop_pad(otherdims_a{:}) = nop;
nop_tail = nop_ext - nop_pad;
N_FM_nop_tail = ones(size(nop_tail));
nop_tail_vec = iso_3D_mat2vec_Fm(nop_tail,N_FM_nop_tail);


%======================================================
% Bounds Y0
%======================================================
% (Y0)_1
Y0_1_body = norm_1_nu_q_3D(AN*F1N,nu,N_Fm,q);
Y0_1_tail = 1/mu_m*norm_1_nu_q_3D(nop_tail_vec,nu,N_FM_nop_tail ,q);
Y0 = sup(Y0_1_body+Y0_1_tail);
%======================================================
% Bounds Z0
%======================================================
AAdagger = AN*DFN;

numel_N = sum(sum(sum(N_Fm)));
B11 = eye(size(AN)) - AAdagger(1:numel_N,1:numel_N);

Z0 = sup(norm_B_nu_q_3D(B11,nu,q,N_Fm));
%======================================================
% Bounds Z1
%======================================================
Id_omega_pad = intval(zeros(size(a_pad))); 
Id_omega_pad(1) = 1;
DNa_ext = N_op_lin_int(Nop_coeff,a_pad,Id_omega_pad);
M = size(DNa_ext)-1;
dim = length(M);
s = arrayfun(@(k) 0:k, M(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if q> 0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^q;
    end
else
    K_sq = intval(ones(size(k{1})));
end
KL_DNa = K_sq.*DNa_ext;


z_1_tilda = reshape(Psi_k_int(DNa_ext ,nu,Na),[],1);
if q == 0
    K = intval(eye(Na_num,Na_num));
else
    K = intval(zeros(Na_num,Na_num));
    for i  = 1:Na_num
        K(i,i) =  sum((iso_coeff(i,Na).*L).^q);
    end
end

Nb = length(a_vect);
K1 = intval(zeros(Nb,numel(mid(a))));
for i = 1:numel(mid(a))
    K1(:,i) = iso_3D_mat2vec_Fm(reshape(K(:,i),Na+1), N_Fm);
end
K2 = intval(zeros(Nb,Nb));
for i = 1:Nb
    K2(i,:) = transpose(iso_3D_mat2vec_Fm(reshape(transpose(K1(i,:)),Na+1), N_Fm));
end

K = K2;






z_1_tilda = iso_3D_mat2vec_Fm(reshape(z_1_tilda,Na+1),N_Fm);

Z1_11_body = norm_1_nu_q_3D(abs(AN)*K*z_1_tilda,nu,N_Fm,q);

N_temp = ones(size(KL_DNa));

KL_DNa_vec =iso_3D_mat2vec_Fm(KL_DNa,N_temp );
Z1_11_tail = 1./mu_m*norm_1_nu_q_3D(KL_DNa_vec,nu,N_temp,q);
Z1 = sup(Z1_11_body + Z1_11_tail);
%======================================================
% Bounds Z2
%======================================================


norm_A11 = max([norm_B_nu_q_3D(AN*K,nu,q,N_Fm),abs(K(end))./mu_m]); 

% ||z_1||
norm_a = norm_1_nu_q_3D(a_vect,nu,N_Fm,q);
norm_z1 = intval(zeros(1,length(Nop_coeff)-2));
for i = 2:(length(Nop_coeff)-1)
    for j = 0:i-2
        for k = 0:j
            norm_z1(end-k) = norm_z1(end-k) +  i*abs(Nop_coeff(end,i+1))*norm_a^(i-j-2)*norm_a^(j-k)*nchoosek(j,k);
        end
    end
end

Z2 = sup(norm_A11*norm_z1);
end

