function mu_ast = find_mu(lambda,m,L)
%
m1 = m(1); m2 = m(2);
N1 = 15; N2 = 15;
L1 = L(1); L2 = L(2);
%
[kx,ky] = ndgrid(0:N1-1,0:N2-1);
%
%
mu_k = lambda - (1-(kx*L1).^2-(ky*L2).^2).^2 ;

mu_k(1:m1,1:m2) = 1.1 * min(mu_k,[],'all');
if all(mu_k(:)<0)
  mu_ast = max(mu_k,[],'all');
else
  mu_ast = NaN;
end