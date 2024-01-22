function [mu_k] = mu_k_SH_sseig(N,lambda,L)
dim = length(N);
s = arrayfun(@(k) 0:k, N(1:end), 'UniformOutput', false); 
k = cell(1, dim);
[k{:}] = ndgrid(s{:});

if exist('intval','file') && isintval(lambda)
    K_sq = intval(zeros(size(k{1})));
else
    K_sq = zeros(size(k{1}));
end

for i = 1:dim
    K_sq = K_sq + (k{i}*L(i)).^2;
end
mu_k = lambda - (1-K_sq).^2;
end

