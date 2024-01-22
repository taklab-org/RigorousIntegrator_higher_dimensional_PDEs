function [F] = F_lin_fin_dim_int2(mu_k,Nop_coeff,a_bar,c,h,e_j)
N = size(c)-1;

for i  = 1:length(N)
    otherdims_test{i} = 1:N(i)+1;
end
a_bar_pad = intval(zeros(N+2));
a_bar_pad(otherdims_test{:}) = a_bar;

c_pad= intval(zeros(N+2));
c_pad(otherdims_test{:}) = c;

mu_k_pad= intval(zeros(N+2));
mu_k_pad(otherdims_test{:}) = mu_k;

e_j_ext =  intval(zeros(N(1:end-1)+2));
e_j_ext(otherdims_test{1:end-1})  = e_j;

% unit_pad = 1*ones(1,length(size(a_bar)));
% a_bar_pad = padarray(a_bar,unit_pad ,'post');
% c_pad = padarray(c,unit_pad ,'post');
% mu_k_pad = padarray(mu_k,unit_pad ,'post');
nop = N_op_lin_int(Nop_coeff,a_bar_pad,c_pad);
% 
% indstring=repmat('1:end-2,',[1 ndims(nop)]);
% indstring(end)=[];
% eval(['nop=nop(' indstring ');']);




psi = h/2*(mu_k_pad.*c_pad + nop);
F = lambda_op_int(c_pad) + T_op_int(psi);

% indstring=repmat('1:end-1,',[1 ndims(nop)]);
% indstring(end)=[];
% eval(['F=F(' indstring ');']);

otherdims = repmat({':'},1,ndims(F));
F(otherdims{:}, 1) =  - e_j_ext;
for j = 0:N(end)
    if j == 0
        F(otherdims{:}, 1) = F(otherdims{:}, 1) + c_pad(otherdims{:}, j+1);
    else
        F(otherdims{:}, 1) = F(otherdims{:}, 1) + 2*((-1)^j)*c_pad(otherdims{:}, j+1);
    end
end

F = F(otherdims_test{:});
end

