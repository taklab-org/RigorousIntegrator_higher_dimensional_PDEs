function [norm_B] = norm_B_nu_q_int_NtoM(B,nu,q,N_Fm,N_Fm_ext)
N_lin_numel = sum(N_Fm+1);
N_col_numel = sum(N_Fm_ext+1);
[N1,N2] = omega12(N_Fm,N_Fm_ext,nu);

if N_col_numel - N_lin_numel >= 0
    i_ind = abs(N_col_numel - N_lin_numel);
    N = N_Fm;
    omega_k = N2;
else
    i_ind = min(N_col_numel, N_lin_numel);
    N_col_new = [N_Fm_ext,-ones(1, length(N_Fm) - length(N_Fm_ext) )];
    N = N_Fm - N_col_new  -1 ;
    omega_k = N1;
end

Norm_colum = intval(zeros(i_ind,1));
for i = 1:i_ind
    B_temp = B(:,i);
    B_temp = iso_vec2mat_Fm_int(B_temp,N);
    Norm_colum(i) = norm_1_nu_q_int(B_temp,nu,q)./(omega_k(i));
end

% omega_k= intval(zeros(1,prod(N+1)));
% for i = 1:prod(N+1)
%     k = iso_coeff(i,N);
%     omega_k(i) = prod(2.^double(k>0))*nu^(sum(k))*(sum(k)+1)^q;
% end
% omega_k = iso_mat2vec_Fm_int(iso_vec2mat_Fm_int(transpose(omega_k),N_lin),N_lin);
% norm_B = max(Norm_colum./omega_k);
norm_B = max(Norm_colum);
end

