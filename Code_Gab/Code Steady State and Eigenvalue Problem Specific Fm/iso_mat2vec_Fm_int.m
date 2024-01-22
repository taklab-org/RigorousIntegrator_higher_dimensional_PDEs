function [b] = iso_mat2vec_Fm_int(a,N_Fm)
N1 = length(N_Fm);
numel_N = sum(N_Fm+1);
b = intval(zeros(numel_N,1));
k = 1;
for i = 1:N1
    for j = 0:N_Fm(i)
        b(k) = a(j+1,i);
        k = k+1;
    end
end
end

