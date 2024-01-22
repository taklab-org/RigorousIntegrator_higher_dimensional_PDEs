function [DF] = DF_steady_state(a,L,mu_k,Nop_coeff,q,N_Fm)
Nb = length(a);
if length(L) == 2
    a = iso_vec2mat_Fm(a,N_Fm);
    mu_k = iso_vec2mat_Fm(mu_k,N_Fm);
elseif length(L) == 3
    a = iso_3D_vec2mat_Fm(a,N_Fm);
    mu_k = iso_3D_vec2mat_Fm(mu_k,N_Fm);
end

N = size(a)-1;
P = size(Nop_coeff,2)-1;
N_ext = P*N;
a_ext = zeros(P*N+1);

for i  = 1:length(N)
    otherdims_test{i} = 1:N(i)+1;
end
a_ext(otherdims_test{:}) = a;


 DF = diag(reshape(mu_k,[],1));
% DF = zeros(numel(a),numel(a));

c = zeros(size(a_ext)); c(1) = 1;
nop = N_op_lin(Nop_coeff,a_ext,c);
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
    K_sq = ones(size(k{1}));
end
nop = nop;


for i = 1:numel(a)
    for j = 1:numel(a)
        k_1 = iso_coeff(i,N);
        k_2 = iso_coeff(j,N);
%         if i == j && i == 9
%             disp('stop here')
%         end

        M = sym_dif(k_2);
        K_q = K_sq(iso_coeff_other(k_1,N_ext));

        for m = 1:size(M,1)
            if K_q ~= 0
                k_out = abs(k_1 - M(m,:) );
                DF(i,j) = DF(i,j) + K_q*nop(iso_coeff_other(k_out,N_ext));
            end
        end
    end
end

DF_temp1 = zeros(Nb,numel(a));
for i = 1:numel(a)
    if length(L) == 2
        DF_temp1(:,i) = iso_mat2vec_Fm(reshape(DF(:,i),N+1), N_Fm);
    elseif length(L) == 3
        DF_temp1(:,i) = iso_3D_mat2vec_Fm(reshape(DF(:,i),N+1), N_Fm);
    end
end
DF_temp2 = zeros(Nb,Nb);
for i = 1:Nb
    if length(L) == 2
        DF_temp2(i,:) = transpose(iso_mat2vec_Fm(reshape(transpose(DF_temp1(i,:)),N+1), N_Fm));
    elseif length(L) == 3
        DF_temp2(i,:) = transpose(iso_3D_mat2vec_Fm(reshape(transpose(DF_temp1(i,:)),N+1), N_Fm));
    end
end

DF = DF_temp2;
end