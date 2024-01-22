function [C,lambda] = Semigroup_estimate_SH2D(a,mu_inf,nu,q,mu_k,Nop_coeff,P,N_Fm,N_Fm_ext,p_r,DF11,N12,N21)
%% Parameters and Variables

a_vect = a;
mu_k_vec = mu_k;
a = iso_vec2mat_Fm_int(a,N_Fm);
mu_k = iso_vec2mat_Fm_int(mu_k,N_Fm);

Na = size(a)-1; % Size of the coefficents 
P_inv = inv(P); % Rigorous inververse on P

C1 = max([norm_B_nu_q_int(P,nu,Na,q,N_Fm),1]);
C2 = max([norm_B_nu_q_int(P_inv,nu,Na,0,N_Fm),1]);

Ma = size(a);
for i = 1:length(Ma)
    otherdims_a{i} = 1:Ma(i); 
end

% Padding of a to match the biggest index of the convolutions
Ma = (size(Nop_coeff,2)-1)*Na; 
a_pad = intval(zeros(Ma+1)); 
a_pad(otherdims_a{:}) = a;

% tic
% [DF11,N12,N21,~] = DF_steady_state_int_ext(a_temp,L,mu_k_temp,Nop_coeff,q,N_Fm_ext,N_Fm);
% toc

A = P_inv*DF11*P; 
eigenvalues = diag(A);
mu_1 = eigenvalues(1);
A_nodiag = A.*-(eye(height(A))-1); % Removing the diagonal term

% delta constant
delta_a = norm_B_nu_q_int(A_nodiag,nu,Na,q,N_Fm) + p_r*C1*C2;

D12 = P_inv*N12;
D21 = N21*P;

delta_b = norm_B_nu_q_int_NtoM(D12,nu,q,N_Fm,N_Fm_ext) + p_r*C2;
delta_c = norm_B_nu_q_int_NtoM(D21,nu,q,N_Fm_ext,N_Fm) + p_r*C1;

[a2_short,a2] = convapbqcos(a,1,a,1);
% delta_b = norm_B_nu_int(P_inv,nu,Na)*3*norm_1_nu_q_int(a2,nu,0);

a2_short = iso_mat2vec_Fm_int(a2_short,N_Fm);


F_m = N_Fm(1)+1;
for i = 1:length(N_Fm)    
    F_m = min([N_Fm(i)+i,F_m]);
end


% delta_b = 3*norm_1_nu_q_int(iso_vec2mat_Fm_int(P_inv*a2_short,N_Fm),nu,q);%/(2*nu^F_m) ;


delta_d = 3*norm_1_nu_q_int(a2,nu,q) + p_r;

% omega_k= intval(zeros(prod(Na+1),1));
% for i = 1:prod(Na+1)
%     k = iso_coeff(i,Na);
%     omega_k(i) = prod(2.^double(k>0))*nu^(sum(k))*(sum(k)+1)^q;
% end
% omega_k = iso_mat2vec_Fm_int(iso_vec2mat_Fm_int(omega_k,N_Fm),N_Fm);
% delta_c = max(norm_Nonlin_inf./transpose(omega_k));
% delta_c = 3/4*norm_1_nu_q_int(a2,nu,0)/(nu^20);

norm_Lambda_inf_inv = abs(1/mu_inf); 
epsilon = intval(0);
for i = 1:length(eigenvalues)
    epsilon = epsilon + 1/(1 -  norm_Lambda_inf_inv*(delta_d + abs(eigenvalues(i))));
end
epsilon = norm_Lambda_inf_inv*epsilon;

if sup(norm_Lambda_inf_inv*(delta_d + max(abs(eigenvalues)))) >= 1
    error('The first condition of the Semigroup estimate is not satisfied')
elseif mu_inf + (delta_d+epsilon*delta_b*delta_c + epsilon^3*delta_b^2*delta_c^2) >= mu_1
    error('The second condition of the Semigroup estimate is not satisfied')
end

% Constant C

C3 = (1+epsilon*delta_b)^2*(1+epsilon*delta_c)^2;

C = C1*C2*C3;
Delta = epsilon*delta_b*delta_c*(max(1+epsilon*delta_c*(1+epsilon*delta_b),epsilon*delta_b*(epsilon^3*delta_b*delta_c)));
lambda = mu_1 + delta_a + Delta;





end

