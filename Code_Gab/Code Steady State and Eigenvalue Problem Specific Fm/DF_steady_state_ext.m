function [N11,N12,N21,N22] = DF_steady_state_ext(a,L,mu_k,Nop_coeff,q,N_Fm_ext,N_Fm)
Nb = length(a);
a = iso_vec2mat_Fm(a,N_Fm_ext);
mu_k = iso_vec2mat_Fm(mu_k,N_Fm_ext);

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
    DF_temp1(:,i) = iso_mat2vec_Fm(reshape(DF(:,i),N+1), N_Fm_ext);
end
DF_temp2 = zeros(Nb,Nb);
for i = 1:Nb
    DF_temp2(i,:) = transpose(iso_mat2vec_Fm(reshape(transpose(DF_temp1(i,:)),N+1), N_Fm_ext));
end
DF = DF_temp2;

N_small = sum(N_Fm+1);
N_big = sum(N_Fm_ext+1);

N11 = [];
N12 = [];
N21 = [];
N22 = [];


N_Fm_pad = [N_Fm,-ones(1,length(N_Fm_ext) - length(N_Fm)) ];

ell = 1;
j_ind = 0;
for j = N_Fm_ext
    for m = 0:j
        k = 1;
        i_ind = 0;
        for i = N_Fm_ext
            for n = 0:i
                if n <= N_Fm_pad(i_ind+1) &&  m <= N_Fm_pad(j_ind+1)
                    N11 = [N11;DF(k,ell)];
                elseif n <= N_Fm_pad(i_ind+1) &&  m > N_Fm_pad(j_ind+1)
                    N12 = [N12;DF(k,ell)];
                elseif n > N_Fm_pad(i_ind+1) &&  m <= N_Fm_pad(j_ind+1)
                    N21 = [N21;DF(k,ell)];
                elseif n > N_Fm_pad(i_ind+1) &&  m > N_Fm_pad(j_ind+1)
                    N22 = [N22;DF(k,ell)];
                end
                k = k+1;
            end
            i_ind =  i_ind  +1 ;
        end
        ell = ell+1;
    end
    j_ind =  j_ind  +1 ;

end

N11 = reshape(N11,N_small,N_small);
N12 = reshape(N12,N_small,N_big-N_small);
N21 = reshape(N21,N_big-N_small,N_small);
N22 = reshape(N22,N_big-N_small,N_big-N_small);















end