function [norm_B] = norm_B_nu(B,nu,N)

N_numel = prod(N+1);
Norm_colum = zeros(1,N_numel);
for i = 1:N_numel
    B_temp = B(:,i);
    B_temp = reshape(B_temp,N+1);
    Norm_colum(i) = norm_1_nu(B_temp,nu);
end
omega_k= zeros(1,N_numel);
for i = 1:N_numel
    k = iso_coeff(i,N);
    omega_k(i) = prod(2.^double(k>0))*nu^(sum(k));
end
norm_B = max(Norm_colum./omega_k);
end
