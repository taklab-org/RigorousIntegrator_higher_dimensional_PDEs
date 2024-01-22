function [omega_k] = omega_weight(nu,N_Fm,q)

if nargin < 3
    q = [];
end
if isempty(q)
    q = 0;
end

a_ones = intval(ones( sum(N_Fm+1) ,1));
a_ones = iso_vec2mat_Fm_int(a_ones,N_Fm);

Ma = size(a_ones)-1;
[m,n] = meshgrid(0:Ma(2),0:Ma(1));

brac_k = (1+m+n).^q;

alpha = 4*intval(ones(size(a_ones)));
alpha(1,:) = 2;
alpha(:,1) = 2;
alpha(1,1) = 1;

omega_k = brac_k .*alpha.*nu.^(m+n);
omega_k = iso_mat2vec_Fm_int(omega_k,N_Fm);



end

