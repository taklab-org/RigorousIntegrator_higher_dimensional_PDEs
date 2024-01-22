function [Y0,Z0,Z1,Z2] = Bounds(a,phi,lambda,L,mu_k,Nop_coeff,nu,q,mu_m)

% Parameters
Na_num = numel(a);
Nphi_num = numel(phi);

Na = size(a)-1;
Nphi = size(phi)-1;

Mphi = size(phi);
Ma = size(a);
for i = 1:length(Ma)
    otherdims_a{i} = 1:Ma(i);
    otherdims_phi{i} = 1:Mphi(i);

end

Ma = (size(Nop_coeff,2)-1)*Na;
a_pad = zeros(Ma+1);
a_pad(otherdims_a{:}) = a;


phi_pad = zeros(Ma+1);
phi_pad(otherdims_phi{:}) = phi;

% Ma = (size(Nop_coeff,2)-1)*Na; % Size of the extension(e.g. Ma = 3*Na for 2D SH)
% unit_pad = (Ma-Na).*ones(1,length(size(a))); %use to pad a
% a_pad = padarray(a,unit_pad ,'post');% Padded a_bar 
% phi_pad = padarray(phi,unit_pad ,'post');% Padded phi_bar 

% Operators
FN_xbar = F_all(a,phi,lambda,L,mu_k,Nop_coeff);
F1N = FN_xbar(1:Na_num);
F2N = FN_xbar(Na_num+1:Nphi_num+Na_num);
F3N = FN_xbar(end);

DFN = DF_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff);

AN = inv(mid(DFN));
AN11 = AN(1:Na_num,1:Na_num);
AN12 = AN(1:Na_num,Na_num+1:Na_num+Nphi_num);
AN13 = AN(1:Na_num, Na_num+Nphi_num+1 );

AN21 = AN(Na_num+1:Na_num+Nphi_num,1:Na_num);
AN22 = AN(Na_num+1:Na_num+Nphi_num,Na_num+1:Na_num+Nphi_num);
AN23 = AN(Na_num+1:Na_num+Nphi_num, Na_num+Nphi_num+1 );

AN31 = AN(Na_num+Nphi_num+1,1:Na_num);
AN32 = AN(Na_num+Nphi_num+1,Na_num+1:Na_num+Nphi_num);
AN33 = AN(Na_num+Nphi_num+1, Na_num+Nphi_num+1 );

nop_ext = N_op_nonlin(Nop_coeff,a_pad);
M_ext = size(nop_ext)-1;
dim = length(M_ext);
s = arrayfun(@(k) 0:k, M_ext(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1) - 1 ) >0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1) - 1 );
    end
else
    K_sq = ones(size(k{1}));
end
nop_ext = K_sq.*nop_ext;

nop= N_op_nonlin(Nop_coeff,a);
M = size(nop)-1;
dim = length(M);
s = arrayfun(@(k) 0:k, M(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1) - 1 ) >0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1) - 1 );
    end
else
    K_sq = ones(size(k{1}));
end
nop = K_sq.*nop;

nop_pad = zeros(Ma+1);
nop_pad(otherdims_a{:}) = nop;


nop_tail = nop_ext - nop_pad;

Dnop = N_op_lin(Nop_coeff,a,phi);
M = size(Dnop)-1;
dim = length(M);
s = arrayfun(@(k) 0:k, M(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1) - 1 ) > 0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1) - 1 );
    end
else
    K_sq = ones(size(k{1}));
end
Dnop = K_sq.*Dnop;

Dnop_ext = N_op_lin(Nop_coeff,a_pad,phi_pad);
M = size(Dnop_ext)-1;
dim = length(M);
s = arrayfun(@(k) 0:k, M(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1) - 1 ) > 0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1) - 1 );
    end
else
    K_sq = ones(size(k{1}));
end
Dnop_ext = K_sq.*Dnop_ext;

Dnop_pad = zeros(Ma+1);
Dnop_pad(otherdims_a{:}) = Dnop;

Dnop_tail = Dnop_ext - Dnop_pad;
%======================================================
% Bounds Y0
%======================================================
% (Y0)_1
Y0_1_body = norm_1_nu(reshape(AN11*F1N,Na+1),nu) + norm_1_nu(reshape(AN12*F2N,Na+1),nu) + norm_1_nu(reshape(AN13*F3N,Na+1),nu);
Y0_1_tail = 1/mu_m*norm_1_nu(nop_tail,nu);
Y0_1 = Y0_1_body+Y0_1_tail;
% (Y0)_2
Y0_2_body = norm_1_nu(reshape(AN21*F1N,Na+1),nu) + norm_1_nu(reshape(AN22*F2N,Na+1),nu) + norm_1_nu(reshape(AN23*F3N,Na+1),nu);
Y0_2_tail = 1/mu_m*norm_1_nu(Dnop_tail,nu);
Y0_2  = Y0_2_body+Y0_2_tail;
% (Y0)_3
Y0_3= abs(AN31*F1N) + abs(AN32*F2N) + abs(AN33*F3N);
Y0 = max([Y0_1,Y0_2,Y0_3]);
%======================================================
% Bounds Z0
%======================================================
AAdagger = AN*DFN;

B11 = eye(size(AN11)) - AAdagger(1:Na_num,1:Na_num);
B12 = -AAdagger(1:Na_num,Na_num+1:Na_num+Nphi_num);
B13 = -AAdagger(1:Na_num, Na_num+Nphi_num+1 );

B21 = -AAdagger(Na_num+1:Na_num+Nphi_num,1:Na_num);
B22 = eye(size(AN22)) -AAdagger(Na_num+1:Na_num+Nphi_num,Na_num+1:Na_num+Nphi_num);
B23 = -AAdagger(Na_num+1:Na_num+Nphi_num, Na_num+Nphi_num+1 );

B31 = -AAdagger(Na_num+Nphi_num+1,1:Na_num);
B32 = -AAdagger(Na_num+Nphi_num+1,Na_num+1:Na_num+Nphi_num);
B33 = eye(size(AN33)) -AAdagger(Na_num+Nphi_num+1, Na_num+Nphi_num+1 );

Z0_1 = norm_B_nu(B11,nu,Na) + norm_B_nu(B12,nu,Na)+ norm_1_nu(reshape(B13,Na+1),nu);
Z0_2 = norm_B_nu(B21,nu,Na) + norm_B_nu(B22,nu,Na)+ norm_1_nu(reshape(B23,Na+1),nu);
Z0_3 = norm_1_nu(reshape(B31',Na+1),nu)+ norm_1_nu(reshape(B32',Na+1),nu) +abs(B33);

Z0 = max([Z0_1,Z0_2,Z0_3]);
%======================================================
% Bounds Z1
%======================================================
Id_omega_pad = zeros(size(a_pad)); 
Id_omega_pad(1) = 1;
DNa_ext = N_op_lin(Nop_coeff,a_pad,Id_omega_pad);
M = size(DNa_ext)-1;
dim = length(M);
s = arrayfun(@(k) 0:k, M(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1) - 1 ) > 0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1) - 1 );
    end
else
    K_sq = ones(size(k{1}));
end
KL_DNa = K_sq.*DNa_ext;

D2Naphi_ext = N_op_lin_twice(Nop_coeff,a_pad,phi_pad);
M = size(D2Naphi_ext)-1;
dim = length(M);
s = arrayfun(@(k) 0:k, M(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1) - 1 ) >0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1) - 1 );
    end
else
    K_sq = ones(size(k{1}));
end
KL_D2Naphi = K_sq.*D2Naphi_ext;

Id_omega = zeros(size(a)); 
Id_omega(1) = 1;

z_1_tilda = reshape(Psi_k(DNa_ext ,nu,Na),[],1);
z_2_tilda = z_1_tilda + reshape(Psi_k(D2Naphi_ext ,nu,Na),[],1);
if q == 0
    K = eye(Na_num,Na_num);
else
    K = zeros(Na_num,Na_num);
    for i  = 1:Na_num
        K(i,i) =  sum((iso_coeff(i,Na).*L).^q);
    end
end

Z1_11_body = norm_1_nu(reshape(abs(AN11)*K*z_1_tilda,Na+1),nu);
Z1_11_tail = 1./mu_m*norm_1_nu(KL_DNa,nu);
Z1_11 = Z1_11_body + Z1_11_tail;

Z1_12 = norm_1_nu(reshape(abs(AN12)*K*z_2_tilda,Na+1),nu);

Z1_21 = norm_1_nu(reshape(abs(AN21)*K*z_1_tilda,Na+1),nu);

Z1_22_body = norm_1_nu(reshape(abs(AN22)*K*z_2_tilda,Na+1),nu);
Z1_22_tail = 1./mu_m*(norm_1_nu(KL_DNa,nu) + norm_1_nu(KL_D2Naphi,nu) ) + abs(lambda)/mu_m;
Z1_22 = Z1_22_body + Z1_22_tail;

Z1_31 = abs(abs(AN31)*K*z_1_tilda);
Z1_32 = abs(abs(AN32)*K*z_2_tilda);

Z1_1 = Z1_11 + Z1_12;
Z1_2 = Z1_21 + Z1_22;
Z1_3 = Z1_31 + Z1_32;

Z1 = max([Z1_1,Z1_2,Z1_3]);
%======================================================
% Bounds Z2
%======================================================


norm_A11 = max([norm_B_nu(AN11*K,nu,Na),abs(K(end))./mu_m]); 
norm_A12 = norm_B_nu(AN12*K,nu,Na); 

norm_A21 = norm_B_nu(AN21*K,nu,Na); 
norm_A22 = max([norm_B_nu(AN22*K,nu,Na),abs(K(end))./mu_m]); 


norm_A31 = norm_1_nu(reshape((AN31*K)',Na+1),nu);
norm_A32 = norm_1_nu(reshape((AN32*K)',Na+1),nu);

% ||z_1||
norm_a = norm_1_nu(a,nu);
norm_phi = norm_1_nu(phi,nu);
norm_z1 = zeros(1,length(Nop_coeff)-2);
norm_z2 = zeros(1,length(Nop_coeff)-2);
for i = 2:(length(Nop_coeff)-1)
    for ell = 0:i-2
        norm_z2(end-ell) = norm_z2(end-ell) + i*(i-1)*abs(Nop_coeff(end,i+1))*norm_a^(i-2-ell);
    end
    for j = 0:i-2
        for k = 0:j
            norm_z1(end-k) = norm_z1(end-k) +  i*abs(Nop_coeff(end,i+1))*norm_a^(i-j-2)*norm_a^(j-k)*nchoosek(j,k);
            norm_z2(end-k) = norm_z2(end-k) +  i*(i-1)*abs(Nop_coeff(end,i+1))*norm_a^(i-j-2)*norm_phi*norm_a^(j-k)*nchoosek(j,k);
        end
    end
end
norm_z2 = norm_z2+norm_z1;
norm_z2(end) = norm_z2(end) + 2;

Z2_1 = norm_A11*norm_z1 + norm_A12*norm_z2;
Z2_2 = norm_A21*norm_z1 + norm_A22*norm_z2;
Z2_3 = norm_A31*norm_z1 + norm_A32*norm_z2;

[~,k_max] = max([Z2_1(end),Z2_2(end),Z2_3(end)]);

Z2_temp = [Z2_1;Z2_2;Z2_3];
Z2 = Z2_temp(k_max,:);

end

