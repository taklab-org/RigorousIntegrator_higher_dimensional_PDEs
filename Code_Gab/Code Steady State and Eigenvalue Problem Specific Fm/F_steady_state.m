function [F] = F_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm)
if length(L) == 2
    a = iso_vec2mat_Fm(a,N_Fm);
    mu_k = iso_vec2mat_Fm(mu_k,N_Fm);
elseif length(L) == 3
    a = iso_3D_vec2mat_Fm(a,N_Fm);
    mu_k = iso_3D_vec2mat_Fm(mu_k,N_Fm);
end

nop = N_op_nonlin(Nop_coeff,a);
M = size(nop)-1;
dim = length(M);
s = arrayfun(@(k) 0:k, M(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if q > 0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^q;
    end
else
    K_sq = ones(size(k{1}));
end
nop = K_sq.*nop;
F = mu_k.*a+nop;

if length(L) == 2
    F = iso_mat2vec_Fm(F,N_Fm);
elseif length(L) == 3
    F = iso_3D_mat2vec_Fm(F,N_Fm);
end

end
