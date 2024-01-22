function [b] = iso_vec2mat_Fm(a,N_Fm)
N2 = length(N_Fm);
N1 = N_Fm(1)+1;
b = zeros(N1,N2);
k = 1;
for i = 1:N2
    for j = 0:N_Fm(i)
        b(j+1,i) = a(k);
        k = k+1;
    end
end
end

