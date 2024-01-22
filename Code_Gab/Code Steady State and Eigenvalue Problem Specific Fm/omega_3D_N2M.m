function [omega_N,omega_inf] = omega_3D_N2M(N_Fm,N_Fm_ext,nu)

omega_k = omega_weight_3D(nu,N_Fm_ext,0);

omega_N = [];
omega_inf = [];

Nc = size(N_Fm);
N1 = length(Nc);
N_all = sum(N_Fm_ext);
for i =1:N1
    N_all = sum(N_all);
end
 
N_Fm_pad = zeros(size(N_Fm_ext));
Size_N_Fm = size(N_Fm);
N_Fm_pad(1:Size_N_Fm(1),1:Size_N_Fm(2),1:Size_N_Fm(3)) = N_Fm;

indicator = iso_3D_mat2vec_Fm(N_Fm_ext - N_Fm_pad,N_Fm_ext);



for i = 1:N_all
    if indicator(i) == 0
        omega_N = [omega_N;omega_k(i)];
    else
        omega_inf = [omega_inf;omega_k(i)];

    end
end

end
