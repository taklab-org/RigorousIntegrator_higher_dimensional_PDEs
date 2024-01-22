function [brac_n,lambda_hat_n] = constant_eta(lambda_hat,t,q)
% Assuming that the eigenvalues are given into a square matrix for for DC2D
N = size(lambda_hat)-1;
[k2,k1] = meshgrid(0:N(1),0:N(2));

brac_k = (1 + k1 + k2).^q;

RHS = brac_k.*exp(-lambda_hat*t);
[row, col] = find(ismember(RHS, max(RHS(:))));

brac_n = brac_k(row,col);
lambda_hat_n = lambda_hat(row,col);

end

