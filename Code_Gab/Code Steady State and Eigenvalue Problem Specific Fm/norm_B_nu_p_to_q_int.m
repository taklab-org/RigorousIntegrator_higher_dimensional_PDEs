function [norm_B] = norm_B_nu_p_to_q_int(B,nu,N,p,q)

N_numel = prod(N+1);
Norm_colum = intval(zeros(1,N_numel));
for i = 1:N_numel
    B_temp = B(:,i);
    B_temp = reshape(B_temp,N+1);
    Norm_colum(i) = norm_1_nu_q_int(B_temp,nu,q);
end
omega_k= intval(zeros(1,N_numel));
for i = 1:N_numel
    k = iso_coeff(i,N);
    omega_k(i) = prod(2.^double(k>0))*nu^(sum(k))*(sum(k)+1)^p;
end
norm_B = max(Norm_colum./omega_k);
end

