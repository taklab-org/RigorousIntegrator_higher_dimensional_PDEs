function [norm_B] = norm_B_3D_N2M(B,nu,N_Fm,M_Fm)
Nc = size(N_Fm);
N1 = length(Nc);
numel_N = sum(N_Fm);
numel_M = sum(M_Fm);
for i =1:N1
    numel_N = sum(numel_N);
    numel_M = sum(numel_M);
end

if numel_M > numel_N
    numel_delta = numel_M - numel_N;
    [~,omega_inf] = omega_3D_N2M(N_Fm,M_Fm,nu);

    Norm_colum = intval(zeros(numel_delta,1));
    for i = 1:numel_delta
        B_temp = B(:,i);
        Norm_colum(i) = norm_1_nu_q_3D(B_temp,nu,N_Fm,0);
    end
    norm_B = max(Norm_colum./omega_inf);
else
    [omega_fin,omega_inf] = omega_3D_N2M(M_Fm,N_Fm,nu);
    Norm_colum = intval(zeros(numel_M,1));
    for i = 1:numel_M
        Norm_colum(i) = sum(abs(B(:,i) ).*omega_inf);
    end
    norm_B = max(Norm_colum./omega_fin);

end

end
