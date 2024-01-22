function [norm_B] = norm_B_X2Y(B,nu1,nu2,q1,q2,N_Fm)

N_numel = sum(N_Fm+1);
Norm_colum = intval(zeros(N_numel,1));
for i = 1:N_numel
    B_temp = B(:,i);
    B_temp = iso_vec2mat_Fm_int(B_temp,N_Fm);
    Norm_colum(i) = norm_1_nu_q_int(B_temp,nu2,q2);
end
omega_k= intval(omega_weight(nu1,N_Fm,q1));
norm_B = max(Norm_colum./omega_k);
end


