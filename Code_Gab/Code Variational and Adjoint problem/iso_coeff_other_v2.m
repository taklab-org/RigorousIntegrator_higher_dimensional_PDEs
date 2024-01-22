function [k] = iso_coeff_other_v2(N)
k = 1:numnel(N)
dim =length(N);
k = 1;
for i = 1:dim
    if i == 1
         k = k + M(i);
    else
        k = k + prod(N(1:i-1)+1)*M(i);
    end
end
end