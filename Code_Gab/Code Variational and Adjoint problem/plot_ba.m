function [c,x,y,t] = plot_ba(a,t,x,y,L)

[x,y,t] = meshgrid(x,y,t);
a = permute(a,[2 3 1]);
N = size(a);
T = Cheb_Taylor_coeff(N(end)-1);
c = zeros(size(t));
for ell = 1:N(1)
    T_ell = 0;
    for j = 0:ell-1
        T_ell = T_ell + T(ell,j+1).*t.^j;
    end
    for k1 = 1:N(1)
        for k2 = 1:N(2)
            if k1 > 1 && k2> 1
                alpha = 4;
            elseif k1 > 1 && k2 == 1
                alpha = 2;
            elseif k1 == 1 && k2 > 1
                alpha = 2;
            elseif k1 == 1 && k2 == 1
                alpha = 1;
            end
            c = c + alpha*a(k1,k2,ell)*cos(k1*L(1)*x).*cos(k2*L(2)*y).*T_ell;

        end
    end
end

end