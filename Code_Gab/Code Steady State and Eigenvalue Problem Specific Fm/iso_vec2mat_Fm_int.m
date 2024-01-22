function [b] = iso_vec2mat_Fm_int(a,N_Fm)
N1 = length(N_Fm);
b = intval(zeros(max([N_Fm,max(N_Fm+1)])));
k = 1;
for i = 1:N1
    for j = 0:N_Fm(i)
        b(j+1,i) = a(k);
        k = k+1;
    end
end
end

