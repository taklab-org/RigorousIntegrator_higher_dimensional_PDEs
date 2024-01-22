function [psi] = Psi_k(a,nu,Nh)
psi = zeros(Nh+1);
Na = size(a)-1;
omega_k = zeros(Na+Nh+1);
for i = 1:prod(Na+Nh+1)
    k = iso_coeff(i,Na+Nh+1);
    omega_k(i) = prod(2.^double(k>0))*nu^(sum(k));
end

for i = 1:prod(Nh+1)
    for j = 1:prod(Na+Nh+1)
        k1 = iso_coeff(i,Na);
        k2 = iso_coeff(j,Na);
        if max(abs(k1-k2)) <= Na(1) && min(abs(k2)) > Nh(1)
            m = iso_coeff_other(abs(k1-k2),Na);
            psi(i) = max([psi(i), a(m)./omega_k(j)]);
        end
        if max(abs(k1+k2)) <= Na(1) && min(abs(k2)) > Nh(1)
            m = iso_coeff_other(abs(k1-k2),Na);
            psi(i) = max([psi(i), a(m)./omega_k(j)]);

        end

    end
end


end

