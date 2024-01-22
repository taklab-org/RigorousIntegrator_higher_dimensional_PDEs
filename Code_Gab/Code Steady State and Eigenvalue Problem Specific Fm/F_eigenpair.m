function [F] = F_eigenpair(a,phi,lambda,L,mu_k,Nop_coeff)
nop = N_op_lin(Nop_coeff,a,phi);
M = size(nop)-1;
dim = length(M);
s = arrayfun(@(k) 0:k, M(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1) - 1 ) > 0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1) - 1 );
    end
else
    K_sq = ones(size(k{1}));
end
nop = K_sq.*nop;
F = (mu_k-lambda).*phi+nop;
F = reshape(F,numel(a),1);
dim = length(size(a));
PC = phi;
for i = 1:dim
    PC = sum(PC);
end
PC = PC-1; % sum(sum(phi)) - 1;sum(sum(phi.^2)) -1;
F = [F;PC];
end

