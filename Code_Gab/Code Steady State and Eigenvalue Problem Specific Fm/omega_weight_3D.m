function [omega_k] = omega_weight_3D(nu,N_Fm,q)

if nargin < 3
    q = [];
end
if isempty(q)
    q = 0;
end

if exist('intval','file') && isintval(nu)
    alpha = 8*intval(N_Fm);
else
    alpha = 8*N_Fm;
end


Ma = size(N_Fm)-1;
[m,n,k] = meshgrid(0:Ma(2),0:Ma(1),0:Ma(3));

brac_k = (1+m+n+k).^q;

alpha(1,:,:) = 4;
alpha(:,1,:) = 4;
alpha(:,:,1) = 4;

alpha(1,1,:) = 2;
alpha(1,:,1) = 2;
alpha(:,1,1) = 2;

alpha(1,1,1) = 1;

omega_k = brac_k .*alpha.*nu.^(m+n+k);
omega_k = iso_3D_mat2vec_Fm(omega_k,N_Fm);



end

