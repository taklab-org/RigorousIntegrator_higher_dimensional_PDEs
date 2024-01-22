function F = F_DC_2D(a,params,N,L) % a is supposed to except zero mode
% input size(a) = (N1*N2-1, 1): vector
% output size(F) = (N1*N2-1, 1): vector form
% 
% params = [lambda,sigma,epsilon];
a0  = params(1);
sigma   = params(2);
epsilon = params(3);
% 
N1 = N(1); N2 = N(2);
L1 = L(1); L2 = L(2);
a_full = [a0;a]; % a_full includes the zero mode
a_full = reshape(a_full,N1,N2);
% 
[kx,ky] = ndgrid(0:N1-1,0:N2-1);

a3 = powerconvcos(a_full,3); % cubic convolution

bkL2 = (kx*L1).^2 + (ky*L2).^2;

mu_k = bkL2 .* (-epsilon^2*bkL2 + 1) - sigma;

F_full = mu_k .* a_full - bkL2 .* a3(1:N1,1:N2);

F_full = reshape(F_full,N1*N2,1); % F_full(1) = F_full(1) - sigma*lambda;

F = F_full(2:end); % remove the zero mode

end