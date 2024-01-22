function [nop_ext] = banach_alg(a,Nop_coeff,q,nop_ext)


N=  size(a)-1;
dim=length(N);

omega = 4;
norm_a = norm_1_nu_q_int(a,omega,q);
omega_vect = intval(zeros(N+1));

for i = 1:numel(mid(a))
    k = iso_coeff(i,N);
    omega_vect(i) = 1./(prod(2.^double(k>0))*omega^(sum(k))*(sum(k)+1)^q);
end
b = midrad(0,sup(omega_vect));


c = intval(zeros(size(a)));
for j=1:dim
    s{j}=1:N(j)+1;
end
for i = 1:(size(Nop_coeff,2)-1)
        if Nop_coeff(end,i+1) ~= 0
                c = c + i*(Nop_coeff(end,i+1)*norm_a^(i-1))*b;
        end
end


nop_ext = intersect(nop_ext,c);
    
end

