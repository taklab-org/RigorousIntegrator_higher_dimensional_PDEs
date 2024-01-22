function [c,k] = newton_F_int(mu_k,Nop_coeff,a_bar,c,h,e_j)
tol = 5e-11;
size_c = size(c);
F = F_lin_fin_dim_int(mu_k,Nop_coeff,a_bar,c,h,e_j);
numel_c = numel(inf(c));
F = reshape(F,[numel_c,1]);
disp(['Before Newton, ||F(c)|| = ',num2str(norm(inf(F),1))])
k_max = 100;
k=0;
while k<=k_max && norm(inf(F),1)> tol
    DF = DF_lin_fin_dim_int(mu_k,Nop_coeff,a_bar,h);
    DF = reshape(DF,[numel_c,numel_c]);
    c = reshape(c,[numel_c,1]);
    c = c - DF\F;
    c = reshape(c,size_c );
    F = F_lin_fin_dim_int(mu_k,Nop_coeff,a_bar,c,h,e_j);
    F = reshape(F,[numel_c,1]);
    disp(['||F(c)|| = ',num2str(norm(inf(F),1))])
    k=k+1;
end
if norm(F)> tol
    error('Did not converged')
end
if isnan(norm(F,1))
    error('Did not converged')
end
display(['||F(a)|| = ',num2str(norm(inf(F),1)),', # Newton iterations = ',num2str(k)])

end

