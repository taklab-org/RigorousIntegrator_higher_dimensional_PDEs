function [Y_0] = Bound_Y0_alone(mu_k,Nop_coeff,a_bar,c,h,e_j,nu,DF_nm,L,lambda)
if ~exist('lambda','var')
    lambda = [];
end
N = size(c)-1; % Size of the dimensions
n = N(end)+1; % Using the same notation as in the paper here, so may need to change some stuff in previous section of the code.
otherdims = repmat({':'},1,length(N)-1); % Used when only working with the Chebyshev indices.
M = (length(Nop_coeff)-1)*(n-1)+1; % Size of the extension(e.g. M = 3n-2 for 2D SH)
unit_pad = zeros(1,length(size(a_bar))); %use to pad a_bar and c.
unit_pad(end) = (M-n+1);
% a_bar_pad = padarray(a_bar,unit_pad ,'post');% Padded a_bar 

Mc = size(c);
Ma = size(a_bar);
for i = 1:length(Ma)
    otherdims_a{i} = 1:Ma(i);
    otherdims_c{i} = 1:Mc(i);
end

c_ext = intval(zeros(Ma));
c_ext(otherdims_c{:}) = c;

otherdims_ext = repmat({':'},1,length(N)-1);

a_bar_pad = intval(zeros([N(1:end-1),M+1]));
a_bar_pad(otherdims_a{:}) = a_bar;

c_pad =  intval(zeros([N(1:end-1),M+1]));
c_pad(otherdims_a{:}) = c_ext;

% c_pad = padarray(c,unit_pad ,'post');% Padded c
%==========================================================================
% Operators
%==========================================================================
numel_c = numel(mid(c));
numel_a = numel(mid(a_bar));
F_nm = reshape(F_lin_fin_dim_int(mu_k,Nop_coeff,a_bar,c,h,e_j,L,lambda),[numel_c,1]); % Truncated F
A_nm =  inv(mid(DF_nm)); %Numerical inverse of DF^(n,m)
T_N_c_bar = T_op_int(h/2*N_op_lin_int(Nop_coeff,a_bar_pad ,c_pad )); %T(N(c_bar))
%==========================================================================
%Bounds Y0 
%==========================================================================
% Bounds over the finite term||A^(n,m)F^(n,m)||
Y_0_1 = norm_1_nu_int(reshape(A_nm*F_nm,size(c)),nu);
% Linear term shifted by Chebyshev
lambda_k_c_k = mu_k.*c;
Y_0_2 = abs(nu^n/(n)*lambda_k_c_k(otherdims{:}, n));
for i = 1:length(N)-1
    Y_0_2 = sum(Y_0_2);
end
% Finite Tail
Y_0_3 = intval(zeros(size(T_N_c_bar)));
for i = n:M
    Y_0_3(otherdims{:}, i+1) = 1/(2*i)*T_N_c_bar(otherdims{:}, i+1) ;
end
Y_0_3 = norm_1_nu_int(Y_0_3,nu);
Y_0 = sup(Y_0_1  + Y_0_2 +  Y_0_3);
end

