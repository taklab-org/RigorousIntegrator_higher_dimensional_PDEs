function [D11,D12,D21] = DF_split_3D(a,L,mu_k,Nop_coeff,q,N_Fm_ext,N_Fm)
Nb = length(a);
if length(L) == 2
a = iso_vec2mat_Fm_int(a,N_Fm_ext);
mu_k = iso_vec2mat_Fm_int(mu_k,N_Fm_ext);
elseif length(L) == 3
    a = iso_3D_vec2mat_Fm(a,N_Fm_ext);
    mu_k = iso_3D_vec2mat_Fm(mu_k,N_Fm_ext);
end

N = size(a)-1;
P = size(Nop_coeff,2)-1;
N_ext = P*N;
a_ext = intval(zeros(P*N+1));

for i  = 1:length(N)
    otherdims_test{i} = 1:N(i)+1;
end
a_ext(otherdims_test{:}) = a;


DF = diag(reshape(mu_k,[],1));
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
    for j = 1:numel(mid(a))
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

DF_temp1 = intval(zeros(Nb,numel(mid(a))));
for i = 1:numel(mid(a))
    if length(L) == 2
        DF_temp1(:,i) = iso_mat2vec_Fm_int(reshape(DF(:,i),N+1), N_Fm_ext);
    elseif length(L) == 3
        DF_temp1(:,i) = iso_3D_mat2vec_Fm(reshape(DF(:,i),N+1), N_Fm_ext);
    end
end
DF_temp2 = intval(zeros(Nb,Nb));
for i = 1:Nb
    if length(L) == 2
        DF_temp2(i,:) = transpose(iso_mat2vec_Fm_int(reshape(transpose(DF_temp1(i,:)),N+1), N_Fm_ext));
    elseif length(L) == 3
        DF_temp2(i,:) = transpose(iso_3D_mat2vec_Fm(reshape(transpose(DF_temp1(i,:)),N+1), N_Fm_ext));
    end
end

DF = DF_temp2;

% Splitting the derivatve
Nc = size(N_Fm);
N1 = length(Nc);
N_small = sum(N_Fm);
for i =1:N1
    N_small = sum(N_small);
end
 
N_big = Nb - N_small;

D11_temp = [];
D12_temp = [];
D21_temp = [];
%D22_temp = [];

N_Fm_pad = zeros(size(N_Fm_ext));
Size_N_Fm = size(N_Fm);
N_Fm_pad(1:Size_N_Fm(1),1:Size_N_Fm(2),1:Size_N_Fm(3)) = N_Fm;

indicator = iso_3D_mat2vec_Fm(N_Fm_ext - N_Fm_pad,N_Fm_ext);



for i = 1:Nb
    if indicator(i) == 0
        D11_temp = [D11_temp,DF(:,i)];
    else
        D12_temp = [D12_temp,DF(:,i)];
        D21_temp  = [D21_temp ;DF(i,:)];
        %D22_temp  = [D22_temp ;DF(i,:)];
    end
end
D11 = [];
D12 = [];
D21 = [];
%D22 = [];
for i = 1:Nb
    if indicator(i) == 0
        D11 = [D11;D11_temp(i,:)];
        D21 = [D21,D21_temp(:,i)];
        D12 = [D12;D12_temp(i,:)];
    %else
        %D22 = [D22,D22_temp(:,i)];
    end
end


end