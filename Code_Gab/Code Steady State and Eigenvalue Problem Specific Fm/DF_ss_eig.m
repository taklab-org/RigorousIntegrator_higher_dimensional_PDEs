function [DF] = DF_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff)
DF_eig = DF_eigenpair(a,phi,lambda,L,mu_k,Nop_coeff);
DF_ss = DF_steady_state(a,L,mu_k,Nop_coeff);
D_phiF_ss = zeros([numel(a),numel(phi)]); 
D_lambdaF_ss = zeros(numel(a),1);

N = size(a)-1;
P = size(Nop_coeff,2)-1;
N_ext = P*N;
a_ext = zeros(P*N+1);
phi_ext = zeros(P*N+1);
for i  = 1:length(N)
    otherdims_test{i} = 1:N(i)+1;
end
a_ext(otherdims_test{:}) = a;
phi_ext(otherdims_test{:}) = phi;

D_aF_eig = zeros(numel(phi),numel(a));

c = zeros(size(a_ext)); c(1) = 1;
nop = N_op_lin_twice(Nop_coeff,a_ext,phi_ext);
M = size(nop)-1;
dim = length(M);
s = arrayfun(@(k) 0:k, M(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1) - 1 ) >0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1) - 1 );
    end
else
    K_sq = ones(size(k{1}));
end
% nop = K_sq.*nop;


for i = 1:numel(a)
    for j = 1:numel(a)
        k_1 = iso_coeff(i,N);
        k_2 = iso_coeff(j,N);
%         if k_1 == k_2
%             DF(i,j) = DF(i,j) + mu_k(k_1+1);
%         end

        M = sym_dif(k_2);
        K_q = K_sq(iso_coeff_other(k_1,N_ext));

        for m = 1:size(M,1)
            if K_q ~= 0
                k_out = abs(k_1 - M(m,:) );
                D_aF_eig(i,j) = D_aF_eig(i,j) + K_q*nop(iso_coeff_other(k_out,N_ext));
            end
        end
    end
end
D_aP = zeros(1,numel(a));
DF  = [DF_ss, D_phiF_ss,D_lambdaF_ss; [D_aF_eig;D_aP], DF_eig];
end
