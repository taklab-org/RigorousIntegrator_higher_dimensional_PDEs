function [Nout] = N_op_nonlin(Nop_coeff,a)
Nout = zeros(size(a));
N = size(a)-1;
dim=length(N);
for j=1:dim
    s{j}=1:N(j)+1;
end
for i = 1:(size(Nop_coeff,2)-1)
        if Nop_coeff(end,i+1) ~= 0
                ap = convapbqcos(a,i-1,a,1);
                Nout = Nout + Nop_coeff(end,i+1)*ap;
        end
end
end