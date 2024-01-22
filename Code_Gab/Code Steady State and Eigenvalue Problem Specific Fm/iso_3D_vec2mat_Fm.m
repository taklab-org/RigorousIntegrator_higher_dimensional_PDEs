function [b] = iso_3D_vec2mat_Fm(a,N_Fm_map)
N1 = length(size(N_Fm_map));
numel_N = sum(N_Fm_map);
for i =1:N1
    numel_N = sum(numel_N);
end
if exist('intval','file') && isintval(a(1))
    b = intval(zeros(size(N_Fm_map)));
else
    b = zeros(size(N_Fm_map));
end

k = 1;
for i = 1:numel(N_Fm_map)
    if N_Fm_map(i) == 1
    b(i) = a(k);
    k = k+1;
    end
end

end


