function DF = fin_dif(mu_k,Nop_coeff,a_bar,c,e_j,L,lambda)
if ~exist('lambda','var')
    lambda = [];
end
N = size(c)-1;
h = 10^(-6);
F = reshape(F_lin_fin_dim(mu_k,Nop_coeff,a_bar,c,h,e_j,L,lambda),[numel(c),1]);
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
E = reshape(E,[numel(c),numel(c)]);
DF = zeros(numel(c),numel(c));
for n = 1:numel(c)
    ch = c + reshape(h.*E(:,n),N+1);
    Fh = reshape(F_lin_fin_dim(mu_k,Nop_coeff,a_bar,ch,h,e_j,L,lambda),[numel(c),1]);
    DF(:,n)= (Fh-F)/h;
end
DF = reshape(DF,[N+1,N+1]);
end

