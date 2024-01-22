function F = F_SH_2D(a,lambda,N,L)
% input size(a) = (N1*N2, 1): vector
% output size(F) = (N1*N2, 1): vector form
%
N1 = N(1); N2 = N(2);
L1 = L(1); L2 = L(2);
a = reshape(a,N1,N2); % size(a) = [N1,N2];
%
[kx,ky] = ndgrid(0:N1-1,0:N2-1);
%
a3 = powerconvcos(a,3);
%
mu_k = lambda - (1-(kx*L1).^2-(ky*L2).^2).^2 ;
%
F = mu_k.*a - a3(1:N1,1:N2);
%
F = reshape(F,N1*N2,1);
%
end