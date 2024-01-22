function [G_0,C,rho,epsilon,delta] = Gershgorin_General(a,P,mu_k,Nop_coeff,L,r,nu,mu_m,q)
%% Definitions
% - For simplicity, we denote the matrix D_tilde = P^-1 DF(a_tilde) P
%   by the 4 parts D_tilde = [A, B ; C , E]
% - R_j represents the Gershgorin radius of associate with the jth eig.
% - G_circ_j represents the Gershgorin circle of associate with the jth eig.

%% Parameters and Variables
Na = size(a)-1; % Size of the coefficents 
P_inv = inv(P); % Rigorous inververse on P
norm_a_nu = norm_1_nu_q_int(a,nu,q); % Norm of a 

% Setting the dimensions of a
Ma = size(a);
for i = 1:length(Ma)
    otherdims_a{i} = 1:Ma(i); 
end

% Padding of a to match the biggest index of the convolutions
Ma = (size(Nop_coeff,2)-1)*Na; 
a_pad = intval(zeros(Ma+1)); 
a_pad(otherdims_a{:}) = a;

% Prealocation of some variable
Nap_i_ext = intval(zeros(numel(mid(a_pad)),numel(mid(a))));
R_j_B_body = intval(zeros(size(diag(P))));
norm_v_i = intval(zeros(size(diag(P))));
norm_Nonlin_fin = intval(zeros(size(diag(P))));
norm_Nonlin_inf = intval(zeros(size(diag(P))));

DFa = DF_steady_state_int(a,L,mu_k,Nop_coeff,q); %DF(a_bar)
A = P_inv*DFa*P; 
A_nodiag = A.*-(eye(height(A))-1); % Removing the diagonal term

%% Gershgorin Circles
for i = 1:numel(mid(a))
    % Padded ith Column of P
    v_i = reshape(P(:,i),Na+1);
    v_i_pad = intval(zeros(Ma+1));
    v_i_pad(otherdims_a{:}) = v_i;

    % Padded tuncated nonlinear part of DF(a_bar)*v_i
    Nop = N_op_lin_int(Nop_coeff,a,v_i);
    Nop_pad = intval(zeros(Ma+1));
    Nop_pad(otherdims_a{:}) = Nop;
    Nop_ext = N_op_lin_int(Nop_coeff,a_pad,v_i_pad);
    Nop_ext = banach_alg_vi(a_pad,v_i,Nop_coeff,q,Nop_ext);
    Nop_ext_cut = Nop_ext;
    Nop_ext_cut(otherdims_a{:}) = 0;
    norm_Nonlin_inf(i) = norm_1_nu_q_int(Nop_ext_cut,nu,0);

%     norm_Nonlin_inf(i) = norm_1_nu_q_int(Nop_ext - Nop_pad,nu,0);
    Nap_i_ext(:,i) = reshape(Nop_ext - Nop_pad,[],1);

    % Main part of R_j
    R_j_B_body(i) = norm_1_nu_q_int(Nop_ext - Nop_pad,nu,q);
    norm_v_i(i) = norm_1_nu_q_int(v_i,nu,q); % norm of v_i
    norm_Nonlin_fin(i) = norm_1_nu_q_int(reshape(A_nodiag(:,i),Na+1),nu,0);
end

% Remainder term from a_tilde  = a_bar + r t with t in the unit ball
p_r = intval(0);
R_inf = intval(0);
C_inf_r = intval(0);
norm_epsilon_2 = intval(0);
for n = 2:length(Nop_coeff)-1
    p_r = p_r + abs(Nop_coeff(n+1))*n*((norm_a_nu + r)^(n-1)-norm_a_nu^(n-1));
    R_inf = R_inf + abs(Nop_coeff(n+1))*n*(norm_a_nu^(n-1)); % bound of R_j for j not in J
    C_inf_r = C_inf_r + abs(Nop_coeff(n+1))*n*2*(norm_a_nu + r)^(n-1);
    norm_epsilon_2 = norm_epsilon_2 + abs(Nop_coeff(n+1))*n*(norm_a_nu + r)^(n-1);
end

% Bound over the constant (KL)^q of the Nonlinearity
L_q_max = max(L.^q);
k_star = min(Na);
omega_k_star = L_q_max/(2*nu^k_star);

% Matrix of weight
nu_vect = intval(zeros(numel(mid(a)),1)); 
for i = 1:numel(mid(a))
    k = iso_coeff(i,Ma);
    nu_vect(i) = prod(2.^double(k>0))*nu^(sum(k))*(sum(k)+1)^q;
end

nu_vect_ext = intval(zeros(numel(mid(a_pad)),1)); 
for i = 1:numel(mid(a_pad))
    k = iso_coeff(i,Na);
    nu_vect_ext(i) = prod(2.^double(k>0))*nu^(sum(k))*(sum(k)+1)^q;
end

c2 = max([sum(abs(P_inv)),1]);



% R_j
R_j = transpose(sum(abs([A_nodiag;Nap_i_ext]))) + c2*p_r*(norm_v_i);
R_inf = c2*(p_r + R_inf);

% G_circ_j 
G_circ_j = diag(A) + midrad(0,sup(c2*p_r*(norm_v_i) + R_j)) ; % Center of GC
G_circ_inf = mu_m + midrad(0,sup(C_inf_r + R_inf));

% Verify the proof
if inf(G_circ_j(1)) >  max(sup(G_circ_j(2:end))) && inf(G_circ_j(1)) > sup(G_circ_inf)
    disp("The proof for the Gershgorin circle was succesful")
else
     error('The proof for the Gershgorin circle failed')
end
% Largest eigenvalues and it's interval 
G_0 = G_circ_j(1);

% Constant C
C1 = max([norm_B_nu_q_int(P,nu,Na,q),1]);
C2 = max([norm_B_nu_q_int(P_inv,nu,Na,0),1]);

% norm of epsilon

delta_a =  norm_B_nu_p_to_q_int(A_nodiag,nu,Na,0,q);
delta_b = 0;
delta_c = 0;
delta_d = 3 * norm_a_nu^2;

eps = 0;
C_inf = 1;

norm_Lambda__inf_inv = abs(1/mu_m);

if sup(norm_Lambda__inf_inv*(delta_d + max(abs(G_circ_j )) )) >=1
    error('First condition on the semigroup estimate failed')
end



C3 = (1+eps*delta_b)^2*(1+eps*delta_c)^2*max([1,C_inf]);
Delta = eps*delta_b*delta_c*(1+eps*(2*delta_b+delta_c)+eps^2*delta_b*delta_c*(1+eps*delta_b));

lambda_s = G_0 + C3*delta_a + Delta*max([1,C_inf]);

C = C1*C2*C3;

rho = 0;
epsilon = 0;
delta= 0;

return
