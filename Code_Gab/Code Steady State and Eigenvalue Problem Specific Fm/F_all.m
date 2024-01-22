function [F] = F_all(a,phi,lambda,L,mu_k,Nop_coeff)
F1 = F_steady_state(a,L,mu_k,Nop_coeff);
F2 = F_eigenpair(a,phi,lambda,L,mu_k,Nop_coeff);
F = [F1;F2];
end

