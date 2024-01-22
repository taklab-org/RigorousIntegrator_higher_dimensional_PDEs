function F = F_SH_2D_ext(a,lambda,n,N,L)
% Input: size(a) = (n, N1*N2) : Matrix form
% Output: size(F) = (3*n-2,3*N1-2,3*N2-2) : full tensor
%
N1 = N(1); N2 = N(2);
L1 = L(1); L2 = L(2);
a = reshape(a,n,N1,N2); % size(a) = [n,N1,N2];
%
[kx,ky] = ndgrid(0:3*N1-3,0:3*N2-3);
[~,a3_full] = powerconvcos(a,3);
%
mu_k = lambda - (1-(kx*L1).^2-(ky*L2).^2).^2;
%
if exist('intval','file') && isintval(a(1))
  a_ext = intval(zeros(3*n-2,3*N1-2,3*N2-2));
else
  a_ext = zeros(3*n-2,3*N1-2,3*N2-2);
end
a_ext(1:n,1:N1,1:N2) = a ;
%
F = permute(mu_k.*permute(a_ext,[2,3,1]),[3,1,2]) - a3_full;
%
end