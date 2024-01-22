function [norm_B] = norm_B_nu_q_3D(B,nu,q,N_Fm)

Nc = size(N_Fm);
N1 = length(Nc);
numel_N = sum(N_Fm);
for i =1:N1
    numel_N = sum(numel_N);
end


Norm_colum = intval(zeros(numel_N,1));
for i = 1:numel_N
    B_temp = B(:,i);
    Norm_colum(i) = norm_1_nu_q_3D(B_temp,nu,N_Fm,q);
end
omega_k= omega_weight_3D(nu,N_Fm,q);
norm_B = max(Norm_colum./omega_k);
end

