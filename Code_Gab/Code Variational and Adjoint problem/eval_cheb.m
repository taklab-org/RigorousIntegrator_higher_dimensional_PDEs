function [c] = eval_cheb(a,t)
N = size(a);
c = 0;
T = Cheb_Taylor_coeff(N(end)-1);
for i = 1:N(end)
    Tn = 0;
    for j = 0:(i-1)
        Tn = Tn + T(i,j+1).*t.^j;
    end
    c = c + a(:,:,i).*Tn;
end

end

