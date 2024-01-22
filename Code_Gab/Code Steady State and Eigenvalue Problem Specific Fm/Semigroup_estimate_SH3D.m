function [C,lambda_s] = Semigroup_estimate_SH3D(DF11,N12,N21,P,mu_inf,Norm_DN,nu,N_Fm,N_Fm_ext,p_r)
% Parameters and Variables
P_inv = inv(P);
M_11 = P_inv*DF11*P;
Lambda_N = diag(M_11).*eye(size(M_11));
E11 =  M_11 - diag(diag(M_11));
E12 = P_inv*N12;
E21 = N21*P;
mu_1 = min(abs(diag(Lambda_N)));
C1 = norm_B_nu_q_3D(P,nu,0,N_Fm);
C2 = norm_B_nu_q_3D(P_inv,nu,0,N_Fm);
eigenvalues_N = diag(Lambda_N);

delta_a = norm_B_nu_q_3D(E11,nu,0,N_Fm) + C2*C1*p_r;
delta_b = norm_B_3D_N2M(E12,nu,N_Fm,N_Fm_ext) + p_r*C2;
delta_c = norm_B_3D_N2M(E21,nu,N_Fm_ext,N_Fm) + p_r*C1;
delta_d = Norm_DN;

epsilon = (1/mu_inf)*sum(1./(mu_inf - delta_d-abs(eigenvalues_N)));

Cond_1 = 1/mu_inf*(delta_d + max(abs(eigenvalues_N)));
Cond_2 = -mu_inf + (delta_d + epsilon*delta_b*delta_c*(1+epsilon^2 * delta_b*delta_c ));
if Cond_1 >= 1
    error('The first condition is not satisfied. You need a bigger gap between \lambda_N and lambda_N+1')
elseif Cond_2 >= -mu_1
    error('The second condition is not satisfied. You need a bigger pseudo-diagonalization')
end

C_s = ( 1 + epsilon*delta_b )^2 * (1 + epsilon*delta_c)^2 ;
Delta = epsilon*delta_b*delta_c*max([1+epsilon*delta_c*(1+epsilon*delta_b),epsilon*delta_b*(2+epsilon^2*delta_b*delta_c)]);
lambda_s = mu_1 - C_s*delta_a - Delta;
C = C_s*C1*C2;
end
