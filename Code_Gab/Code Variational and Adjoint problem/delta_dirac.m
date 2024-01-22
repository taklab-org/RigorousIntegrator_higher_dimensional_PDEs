function [v] = delta_dirac(x)
v = ones(size(x));
for i = 1:length(x)
    if x(i) == 0
    v(i) = 0;
    end
end

end

