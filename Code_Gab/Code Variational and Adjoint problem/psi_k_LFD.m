function [psi] = psi_k_LFD(mu_k,Nop_coeff,a_bar,c,h)
nop = N_op_lin(Nop_coeff,a_bar,c);
psi = h/2*(mu_k.*c + nop);
end

