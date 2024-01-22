function [mu_k] = mu_k_dcCH_sseig(N,epsilon,sigma,L)
dim = length(N);
s = arrayfun(@(k) 0:k, N(1:end), 'UniformOutput', false); 
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
K_sq = zeros(size(k{1}));
for i = 1:dim
    K_sq = K_sq + (k{i}*L(i)).^2;
end
mu_k = (K_sq).*(-epsilon^2*K_sq +1)-sigma;
end

