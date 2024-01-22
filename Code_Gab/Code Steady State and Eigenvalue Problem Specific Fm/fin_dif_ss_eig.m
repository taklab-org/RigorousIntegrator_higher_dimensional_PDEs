function DF = fin_dif_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff)
phi_lambda = phi; phi_lambda(1) = lambda;
x = [a,phi_lambda];
size_x = size(x);
size_a = size(a);
size_phi = size(phi);
N = size(x)-1;
h = 10^(-6);
F = reshape(F_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff),[numel(a)+numel(phi),1]);
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
E = reshape(E,[numel(x),numel(x)]);
DF = zeros(numel(x),numel(x));
for n = 1:numel(x)
    xh = x + reshape( h.*E(:,n),size_x);
    ah = xh(:,1:size_a(2));
    phih = xh(:,size_a(2)+1:end);
    lambdah = phih(1);
    phih(1) = phi(1);
    Fh = reshape(F_ss_eig(ah,phih,lambdah,L,mu_k,Nop_coeff),[numel(x),1]);
    DF(:,n)= (Fh-F)/h;
end
DF = reshape(DF,[size_x ,size_x ]);
end

