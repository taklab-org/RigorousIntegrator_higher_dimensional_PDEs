function F = F_SH_3D(a,lambda,N,L)
% input size(a) = [N1*N2*N3, 1] : vector
% output size(F) = [N1*N2*N3, 1] : vector
N1 = N(1); N2 = N(2); N3 = N(3);
L1 = L(1); L2 = L(2); L3 = L(3);
a = reshape(a,N1,N2,N3); % size(a) = [N1,N2,N3];
% 
[kx,ky,kz] = ndgrid(0:N1-1,0:N2-1,0:N3-1);
% 
a3 = powerconvcos(a,3);
% 
mu_k = lambda - (1-(kx*L1).^2-(ky*L2).^2-(kz*L3).^2).^2 ;
% 
F = mu_k.*a - a3(1:N1,1:N2,1:N3);
% 
F = reshape(F,N1*N2*N3,1); 
% 
end