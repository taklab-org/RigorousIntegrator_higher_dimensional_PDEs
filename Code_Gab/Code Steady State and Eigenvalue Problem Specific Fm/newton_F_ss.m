function [a,k] = newton_F_ss(a,L,mu_k,Nop_coeff,q,N_Fm)

if length(L) == 2

    tol = 5e-14;
    size_a = size(a);
    F = F_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
    disp(['Before Newton, ||F(c)|| = ',num2str(norm(F,1))])
    k_max = 100;
    k=0;
    while k<=k_max && norm(F,1)> tol
        DF = DF_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
        a = a - DF\F;
        F = F_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
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


elseif length(L) == 3
    tol = 5e-14;
    F = F_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
    disp(['Before Newton, ||F(c)|| = ',num2str(norm(F,1))])
    k_max = 100;
    k=0;
    while k<=k_max && norm(F,1)> tol
        DF = DF_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
        a = a - DF\F;
        F = F_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
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

