function F = F_SH_3D_ext(a,lambda,n,N,L)
% Input: size(a) = [n, N1*N2*N3] : Matrix form
% Output: size(F) = [3*n-2,3*N1-2,3*N2-2,3*N3-2] : full 4D tensor
% 
N1 = N(1); N2 = N(2); N3 = N(3);
L1 = L(1); L2 = L(2); L3 = L(3);
a = reshape(a,n,N1,N2,N3); % size(a) = [n,N1,N2,N3];
% 
[kx,ky,kz] = ndgrid(0:3*N1-3,0:3*N2-3,0:3*N3-3);
% 
[~,a3_full] = powerconvcos(a,3);
% 
mu_k = lambda - (1-(kx*L1).^2-(ky*L2).^2-(kz*L3).^2).^2 ;
% 
if exist('intval','file') && isintval(a(1))
  a_ext = intval(zeros(3*n-2,3*N1-2,3*N2-2,3*N3-2));
else
  a_ext = zeros(3*n-2,3*N1-2,3*N2-2,3*N3-2);
end
% 
a_ext(1:n,1:N1,1:N2,1:N3) = a ;
% 
F = permute(mu_k.*permute(a_ext,[2,3,4,1]),[4,1,2,3]) - a3_full;
% 
end