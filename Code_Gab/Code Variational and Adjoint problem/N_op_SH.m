function [Nout] = N_op_SH(Nop_coeff,a)
Nout = zeros(size(a));
for i = find(~Nop_coeff)
    Nout = Nout + Nop_coeff(i+1)*powerconvcos(a,i);
end
end

