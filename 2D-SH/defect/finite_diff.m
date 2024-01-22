function DF = finite_diff(a,lambda,N,L)

h = 1e-6;
m = length(a);
E = eye(m);
DF = zeros(m);
for j=1:m
    ah = a+h*E(:,j);
    DF(:,j) = (F_SH_2D(ah,lambda,N,L)-F_SH_2D(a,lambda,N,L))/h;
end

end