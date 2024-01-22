function F = vec_F_SH_2D(a,lambda,m,n)
a = reshape(a,m,n);
F = F_SH_2D(a,lambda);
F = reshape(F,m*n,1);

