function [DF] = Dnop_non_lin(a,L,Nop_coeff,q)

N = size(a)-1;
P = size(Nop_coeff,2)-1;
N_ext = P*N;
a_ext = intval(zeros(P*N+1));

for i  = 1:length(N)
    otherdims_test{i} = 1:N(i)+1;
end
a_ext(otherdims_test{:}) = a;


DF = intval(zeros(numel(mid(a)),1));
% DF = zeros(numel(a),numel(a));

c = intval(zeros(size(a_ext))); c(1) = 1;
nop = N_op_lin_int(Nop_coeff,a_ext,c);
M = size(nop)-1;
dim = length(M);
s = arrayfun(@(k) 0:k, M(1:end), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if q>0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^q;
    end
else
    K_sq = intval(ones(size(k{1})));
end

for i = 1:numel(mid(a))
        k_1 = iso_coeff(i,N);
        k_2 = iso_coeff(i,N);
%         if i == j && i == 9
%             disp('stop here')
%         end

        M = sym_dif(k_2);
        K_q = K_sq(iso_coeff_other(k_1,N_ext));

        for m = 1:size(M,1)
            if K_q ~= 0
                k_out = abs(k_1 - M(m,:) );
                DF(i) = DF(i) + K_q*nop(iso_coeff_other(k_out,N_ext));
            end
        end
end
end