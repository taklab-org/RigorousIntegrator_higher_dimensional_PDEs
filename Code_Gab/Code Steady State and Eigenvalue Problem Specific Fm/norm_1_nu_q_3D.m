function [n] = norm_1_nu_q_3D(a,nu,N_Fm,q)
omega = omega_weight_3D(nu,N_Fm,q);
n = sum(abs(a).*omega);
end

