function DF = fin_dif_ss(a,L,mu_k,Nop_coeff,q,N_Fm)

N = size(a)-1;
h = 10^(-6);
F = F_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm);
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
    ah = a + reshape(h.*E(:,n),N+1);
    Fh = F_steady_state(ah,L,mu_k,Nop_coeff,q,N_Fm);
    DF(:,n)= (Fh-F)/h;
end
end

