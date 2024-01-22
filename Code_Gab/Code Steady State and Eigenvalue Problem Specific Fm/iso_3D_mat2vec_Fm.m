function [b] = iso_3D_mat2vec_Fm(a,N_Fm_map)


Na = size(a);
Nc = size(N_Fm_map);
N1 = length(Nc);
N_max = min([Na;Nc]);
numel_N = sum(N_Fm_map);
for i =1:N1
    numel_N = sum(numel_N);
end

if exist('intval','file') && isintval(a(1))
    a_ext = intval(zeros(Nc));
    b = intval(zeros(numel_N,1));
else
    a_ext = zeros(Nc);
    b = zeros(numel_N,1);
end
a_ext(1:N_max(1),1:N_max(2),1:N_max(3)) = a(1:N_max(1),1:N_max(2),1:N_max(3));


k = 1;
for i = 1:numel(N_Fm_map)
    if N_Fm_map(i) == 1
    b(k) = a_ext(i);
    k = k+1;
    end
end
end

