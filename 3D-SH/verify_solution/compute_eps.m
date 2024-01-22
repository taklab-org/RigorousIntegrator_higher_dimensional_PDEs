function eps=compute_eps(ba0,a0,nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input is the one-sided in Chebyshev %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chebfunpref.setDefaults('factory');
% index = 1;
N = size(ba0);
n = N(1); N1 = N(2); N2 = N(3); N3 = N(4);
ba0matrix = reshape(ba0,n,N1*N2*N3);
ba0_twosided = [ba0matrix(1,:);2*ba0matrix(2:end,:)];
%
ell = 0:n-1;
v = (-1).^ell;
% v_hat = zeros(1,2*N+1);
% v_hat(N+1)=50; v_hat(N)=-25; v_hat(N+2)=-25;
%
eps = wnorm(reshape(v*ba0_twosided,N1,N2,N3)-a0,nu);