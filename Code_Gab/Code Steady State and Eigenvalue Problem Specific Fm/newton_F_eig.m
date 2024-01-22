function [phi,lambda,k] = newton_F_eig(a,phi,lambda,L,mu_k,Nop_coeff)

phi_1 = phi(1);
tol = 5e-12;
size_a = size(a);
F = F_eigenpair(a,phi,lambda,L,mu_k,Nop_coeff);
disp(['Before Newton, ||F(c)|| = ',num2str(norm(F,1))])
k_max = 100;
k=0;
while k<=k_max && norm(F,1)> tol
    DF = DF_eigenpair(a,phi,lambda,L,mu_k,Nop_coeff);
    x = [reshape(phi,[numel(a),1]);lambda];
    x = x - DF\F;
    phi = reshape(x(1:end-1),size_a ); lambda = x(end);
    F = F_eigenpair(a,phi,lambda,L,mu_k,Nop_coeff);
    disp(['||F(c)|| = ',num2str(norm(F,1))])
    k=k+1;
end
if norm(F)> tol
    error('Did not converged')
end
if isnan(norm(F,1))
    error('Did not converged')
end
display(['||F(a)|| = ',num2str(norm(F,1)),', # Newton iterations = ',num2str(k)])

end

