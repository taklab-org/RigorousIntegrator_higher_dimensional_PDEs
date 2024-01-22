function DF = fin_dif_eig(a,phi,lambda,L,mu_k,Nop_coeff)

N = size(phi)-1;
h = 10^(-6);
F = reshape(F_eigenpair(a,phi,lambda,L,mu_k,Nop_coeff),[numel(a),1]);
E = zeros([N+1,N+1]);
if length(N) == 3
    for n = 0:N(1)
        for m = 0:N(2)
            for ell = 0:N(3)
                E(n+1,m+1,ell+1,n+1,m+1,ell+1) = 1;
            end
        end
    end
elseif length(N) == 2
    for n = 0:N(1)
        for m = 0:N(2)
            E(n+1,m+1,n+1,m+1) = 1;
        end
    end
end
E = reshape(E,[numel(a),numel(a)]);
DF = zeros(numel(a),numel(a));
for n = 1:numel(a)
    if n == 1
        lambdah = lambda + h;
        Fh = reshape(F_eigenpair(a,phi,lambdah,L,mu_k,Nop_coeff),[numel(a),1]);
    else
        phih = phi + reshape(h.*E(:,n),N+1);
        Fh = reshape(F_eigenpair(a,phih,lambda,L,mu_k,Nop_coeff),[numel(a),1]);
    end
    DF(:,n)= (Fh-F)/h;
end
DF = reshape(DF,[N+1,N+1]);
end

