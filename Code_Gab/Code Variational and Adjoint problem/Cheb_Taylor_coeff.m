function [T] = Cheb_Taylor_coeff(N)
T = zeros(N+1);
T(1,1) = 1;
if N >= 1
    T(2,2) = 1;
    if N > 1
        for k = 2:N
            for j = 0:N
                if j == 0
                    if mod(k,2) == 0
                        T(k+1,j+1) = (-1)^(k/2);
                    end
                else
                    T(k+1,j+1) = 2*T(k,j)- T(k-1,j+1);
                end
            end
        end
    end
end
end

