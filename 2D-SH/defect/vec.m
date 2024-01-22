function f = vec(a)
[m,n] = size(a);
f = reshape(a,m*n,1);