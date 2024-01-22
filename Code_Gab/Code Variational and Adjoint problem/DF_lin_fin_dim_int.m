function [DF] = DF_lin_fin_dim_int(mu_k,Nop_coeff,a_bar,h,c,L,lambda)
if ~exist('lambda','var')
    lambda = [];
end
N = size(a_bar)-1;
P = size(Nop_coeff,2)-2;
% a_bar_ext = padarray(a_bar,L*N,'post');
Mc = size(c);
Ma = size(a_bar);
mu_k_temp = mu_k;
mu_k = intval(zeros(Ma));
for i  = 1:length(N)
    otherdims_test{i} = 1:N(i)+1;
            otherdims_a{i} = 1:Ma(i);
    otherdims_c{i} = 1:Mc(i);
end
mu_k(otherdims_c{:}) = mu_k_temp;
a_bar_ext = intval(zeros((P+1)*N+1));
a_bar_ext(otherdims_test{:}) = a_bar;
Mc = size(c);

% otherdims_ext = repmat({':'},1,ndims(a_bar));
% a_bar_ext(otherdims_ext{:},(N(end)+1):((L+1)*N(end)+1)) = 0 ;

% numel_a_bar = numel(mid(a_bar));
numel_c = numel(mid(c));
N_ext = size(a_bar_ext)-1;
unit = intval(zeros(size(a_bar_ext)));
unit(1) = 1;
nop = N_op_lin_int(Nop_coeff,a_bar_ext,unit);
M = size(nop)-1;
dim = length(M)-1;
s = arrayfun(@(k) 0:k, M(1:end-1), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
K_sq = intval(zeros(size(k{1})));
for ell = 1:dim
    K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1) - 1 );
end
N_rep = ones(size(M));
N_rep(end) = M(end)+1;
K_sq = repmat(K_sq,N_rep);


nop = reshape(nop,[numel(inf(a_bar_ext)),1]);







DF = intval(zeros(numel_c,numel_c));
% tic
% coeff_matrix = (1:numel(a_bar_ext))';
% coeff_matrix = reshape(coeff_matrix,size(a_bar_ext));
% toc

for i = 1:numel_c
    for j = 1:numel_c
        k_1 = iso_coeff(i,Mc-1);
        k_2 = iso_coeff(j,Mc-1);
        if k_1(end) == 0 && isequal(k_1(1:end-1),k_2(1:end-1))
            if k_2(end)  == 0
                DF(i,j) = DF(i,j) +  1;
            else
                DF(i,j) = DF(i,j) +  2*(-1)^k_2(end);
            end
        elseif k_1(end) > 0 && isequal(k_1,k_2)
            DF(i,j) = DF(i,j) + 2*k_1(end);
        elseif k_1(end) > 0 && isequal(k_1(end)-1,k_2(end)) && isequal(k_1(1:end-1),k_2(1:end-1))
            if isempty(lambda) || ~isequal(k_1(1:(end-1)), zeros(size(k_1(1:(end-1)))) )
                DF(i,j) = DF(i,j) - h/2*mu_k(iso_coeff_other(k_1,N));
            end
        elseif k_1(end) > 0 && isequal(k_1(end)+1,k_2(end)) && isequal(k_1(1:end-1),k_2(1:end-1))
            if isempty(lambda) || ~isequal(k_1(1:(end-1)), zeros(size(k_1(1:(end-1)))) )
                DF(i,j) = DF(i,j) + h/2*mu_k(iso_coeff_other(k_1,N));
            end
        end
        if k_1(end) > 0
            if isempty(lambda) || ~isequal(k_1(1:(end-1)), zeros(size(k_1(1:(end-1)))) )
                M = sym_dif(k_2);
                K_q = K_sq(iso_coeff_other(k_1,N_ext));
                for m = 1:size(M,1)
                    K_shift = k_1-M(m,:);

                    K_shift_minus = K_shift;
                    K_shift_minus(end) = K_shift_minus(end)-1;
                    %K_shift_minus = num2cell(abs(K_shift_minus)+1);
                    K_shift_minus = abs(K_shift_minus);

                    K_shift_plus = K_shift;
                    K_shift_plus(end) = K_shift_plus(end)+1;
                    % K_shift_plus = num2cell(abs(K_shift_plus)+1);
                    K_shift_plus = abs(K_shift_plus);

                    DF(i,j) = DF(i,j) + h/2*K_q*nop(iso_coeff_other(K_shift_plus,N_ext)) - ...
                        h/2*K_q*nop(iso_coeff_other(K_shift_minus,N_ext));
                    %                 DF(i,j) = DF(i,j) + h/2*nop(coeff_matrix(K_shift_plus{:})) - ...
                    %                     h/2*nop(coeff_matrix(K_shift_minus{:}));
                end
            end
        end
    end
end
DF = reshape(DF,[Mc,Mc]);

end

