function [n] = norm_1(A)
L = length(size(A));
n = abs(A);
for k = 1:L
    n = max(n);
end
end

