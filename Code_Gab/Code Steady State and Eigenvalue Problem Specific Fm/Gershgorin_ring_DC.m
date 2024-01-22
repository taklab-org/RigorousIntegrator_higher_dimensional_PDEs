function [G_0,C,rho,epsilon,delta] = Gershgorin_ring_DC(a,P,mu_k,Nop_coeff,L,r,nu,mu_m,q)

Na = size(a)-1;
P_inv = inv(P);

Ma = size(a);
for i = 1:length(Ma)
    otherdims_a{i} = 1:Ma(i);
end

Ma = (size(Nop_coeff,2)-1)*Na;
a_pad = intval(zeros(Ma+1));
a_pad(otherdims_a{:}) = a;

Nap_i_ext = intval(zeros(numel(mid(a_pad)),numel(mid(a))));
DFa = DF_steady_state_int(a,L,mu_k,Nop_coeff,q);

A = P_inv*DFa*P;
B = intval(diag(A));
A = A.*-(eye(height(A))-1);

norm_a_nu = norm_1_nu_q_int(a,nu,q);
norm_p_i = intval(zeros(1,numel(mid(a))));
norm_p_i_1 = intval(zeros(1,numel(mid(a))));

a_ext = intval(zeros(3*Na+1));
a_ext(1:Na(1)+1,1:Na(2)+1) = a;

omega_m = 2*nu^(Na(1));




h = intval('0.25')*intval(ones(size(a)));
h(1,:) = intval('0.5');
h(:,1) = intval('0.5');
h(1,1) = intval('1');
h = reshape(h,[],1);

max_a2v = intval(zeros(1,numel(mid(a))));
a2v = intval(zeros(size(a_ext)));
for i = 1:numel(mid(a))
    p_i = reshape(P(:,i),Na+1);
    p_i_pad = intval(zeros(Ma+1));
    p_i_pad(otherdims_a{:}) = p_i;

    Nop = N_op_lin_int(Nop_coeff,a,reshape(P(:,i),Na+1));
    Nop_pad = intval(zeros(Ma+1));
    Nop_pad(otherdims_a{:}) = Nop;

    Nop_ext = N_op_lin_int(Nop_coeff,a_pad,p_i_pad);

    Nap_i_ext(:,i) = reshape(Nop_ext - Nop_pad,[],1);
    norm_p_i(i) = norm_1_nu_q_int(p_i,nu,q);
    norm_p_i_1(i) = norm_1_nu_q_int(p_i,1,q);

    v_i_ext = intval(zeros(3*Na+1));
    v_i_ext(1:Na(1)+1,1:Na(2)+1) = p_i;
    a2v_i = convapbqcos(a_ext,2,v_i_ext,1);
    a2v_i(1:Na(1)+1,1:Na(2)+1) = 0;

    max_a2v(i) = norm_1_nu_q_int(a2v_i,nu,q);
    a2v = a2v + a2v_i*h(i);

end
Epsilon21 = 3*norm_1_nu_q_int(a2v,1,q);

C_2 = max([norm_B_nu_q_int(P_inv,1,Na,0),1]);

R = transpose(sum(abs([A;Nap_i_ext]))+ C_2*(6*r*norm_a_nu+3*r^2)*norm_p_i);
B = B + norm_B_nu_q_int(P_inv,nu,Na,q)*(6*r*norm_a_nu+3*r^2)*transpose(norm_p_i);


B_tail = mu_m + (6*norm_B_nu_int(P,nu,Na)+3)*(norm_a_nu^2 + 2*r*norm_a_nu+r^2);

R_int = midrad(0,sup(R));
G_circ = B + R_int; 

if inf(G_circ(1)) >  max(sup(G_circ(2:end))) && inf(G_circ(1)) > sup(B_tail)
    disp("The proof worked!")
else
    error('The proof failed')
end

G_0 = G_circ(1);
C_1 = max([norm_B_nu_int(P,1,Na),1]);

C = sup(C_1*C_2);

% ||\epsilon||
% epsilon 11

p_r = 2*r*norm_a_nu + r^2;
max_P_inv = transpose(max(transpose(abs(P_inv))));
K = 3*max_P_inv.*norm_p_i;

A_int = inv(P)*intval(DFa)*intval(P);
D_int = diag(A_int);
A_int = A_int.*-(eye(height(A_int))-1);

Ah = reshape(A_int*h,size(a));
Epsilon11 = norm_1_nu_int(Ah + p_r*(reshape(K*abs(h),Na+1)) ,1) ;

% epsilon 12
[~,a2] = convapbqcos(a,1,a,1);


psi_a2 = Psi_k_int(a2,1,Na);
sum_psi_a2_cut = sum(sum(psi_a2(1:Na(1)+1,1:Na(2)+1)));

Epsilon12 = 3*norm_1_nu_int(reshape(max_P_inv,Na+1),1)*(sum_psi_a2_cut+p_r);

% epsilon 21

Epsilon21 = Epsilon21 + 3*p_r*max(norm_p_i);


% epsilon 22
norm_a2_nu = norm_1_nu_int(a2,nu);
Epsilon22 = 3/omega_m*(norm_a2_nu+2*r*norm_a_nu+r^2);
gamma  = sqrt(Epsilon12/Epsilon21);


Epsilon = Epsilon11 + gamma*Epsilon21 + 1/gamma*Epsilon12 + Epsilon22;

% Lambda
lambda_hat = max(D_int);
lambda = abs(sup(lambda_hat + Epsilon)); 

epsilon = 10^-10;
delta = inf(intval(lambda) - intval(epsilon));

p_gamma = [1, mid(3*(norm_a_nu+r)), -delta/C];
rho = max(roots(p_gamma));

end

