function F = F_DC_2D_ext(a,params,n,N,L) % a is supposed to include zero mode
% Input: size(a) = (n, N1*N2) : Matrix form
% Output: size(F) = (3*n-2,3*N1-2,3*N2-2) : full tensor
% 
% params = [lambda,sigma,epsilon];
% lambda  = params(1);
sigma   = params(2);
epsilon = params(3);
% 
N1 = N(1); N2 = N(2);
L1 = L(1); L2 = L(2);
% 
a = reshape(a,n,N1,N2);
% 
[kx,ky] = ndgrid(0:3*N1-3,0:3*N2-3);
[~,a3_full] = powerconvcos(a,3); % cubic convolution
% 
bkL2 = (kx*L1).^2 + (ky*L2).^2;
mu_k = bkL2 .* (-epsilon^2*bkL2 + 1) - sigma;
% 
if exist('intval','file') && isintval(a(1))
  a_ext = intval(zeros(3*n-2,3*N1-2,3*N2-2));
else
  a_ext = zeros(3*n-2,3*N1-2,3*N2-2);
end
a_ext(1:n,1:N1,1:N2) = a;
% 
F = permute(mu_k.*permute(a_ext,[2,3,1]),[3,1,2]) - permute(bkL2.*permute(a3_full,[2,3,1]),[3,1,2]);
% 
end