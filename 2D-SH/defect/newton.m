function a = newton(a,lambda,N,L)

tol=1e-14;
k=0;
F = F_SH_2D(a,lambda,N,L);
disp(['||F|| = ', num2str(norm(F))])
while norm(F)>tol && k<100lambda
    DF = finite_diff(a,lambda,N,L);
    a = a - DF\F_SH_2D(a,lambda,N,L);
    F = F_SH_2D(a,lambda,N,L);
    k=k+1;
    %display(['||F|| = ', num2str(norm(F))])
end

display(['||F|| = ', num2str(norm(F)),', ||DF^(-1)|| = ', num2str(norm(inv(DF)))])

end