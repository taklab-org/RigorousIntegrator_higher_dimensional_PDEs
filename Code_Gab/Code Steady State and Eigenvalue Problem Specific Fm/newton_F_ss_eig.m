function [a,phi,lambda,k] = newton_F_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff)
tol = 5e-12;
F = F_all(a,phi,lambda,L,mu_k,Nop_coeff);
disp(['Before Newton, ||F(c)|| = ',num2str(norm(F,1))])
k_max = 100;
k=0;
while k<=k_max && norm(F,1)> tol
    DF = DF_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff);
    x = [reshape(a,numel(a),1);reshape(phi,numel(phi),1);lambda];
    x = x - DF\F;
    a = reshape(x(1:numel(a)),size(a));
    phi = reshape(x(numel(a)+1:end-1),size(phi));
    lambda = x(end);
    F = F_all(a,phi,lambda,L,mu_k,Nop_coeff);
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

