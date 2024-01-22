function [psi] = psi_z_int(a,nu)
N = size(a)-1;
n = N(end)+1;
numel_a = numel(mid(a));
psi = intval(zeros(numel_a ,1));
b = reshape(abs(a),[numel_a ,1]);
for i = 2:numel_a 
    for j = 1:numel_a 
        k_ell = iso_coeff(i,N);
        k = k_ell(1:end-1);
        ell = k_ell(end);
        k2_ell2 = iso_coeff(j,N);
        ell2 = k2_ell2(end)+ell;
        j2 = k2_ell2(1:end-1);
        if max(abs(k-j2)) <= N(1) || max(abs(k+j2)) <= N(1)
            if ell2 >= n && ell2 <= ell+n-1
                psi(j) = max(psi(j),b(j)/(nu^ell2));
            end
        end
    end
    psi = reshape(psi,size(a));
end


