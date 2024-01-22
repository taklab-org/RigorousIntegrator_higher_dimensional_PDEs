function [r,Y0,Z0,Z1] = rigorous_bounds_lin_SH(c,a_bar,h,e_j,lambda,nu)

% Description:
%--------------------------------------------------------------------------
% - Rigourous computation of the linearization about a_bar for the 
%   Swift-Hohenberg equation:
%   
%             u_t = \lambda u - (1+\Delta)^2 u -u^3.    (1)
%
% - Input:
%         c: numerical approximation of the solution of the linearization.
%     a_bar: approximate solution of (1)
%         h: time interval [0,h]
%       e_j: Canonical basis vector
%    lambda: Parameter of the S-H equation.
%        nu: Weight of the Chebyshev expension.
%
% - Output: 
%          r: Smallest radius of existence of the ball centered at c
%         Y0: Rigorous bound over ||AF(c)||_{\ell_\nu^1}
%         Z0: Rigorous bound over ||Id - AA^\dagger||_{B(X)}
%         Z1: Rigorous bound over ||A(Df(c) - A^\dagger)||_{B(X)}
%--------------------------------------------------------------------------

% Order of coefficients dimensions: 
%--------------------------------------------------------------------------
% To use this code, the dimension of c and a_bar corresponding to the
% Chebyshev expension must be in the last position. If they are in the 
% first position, uncomment the following line to permute the order
% of the dimensions:
%
% L = length(size(c));
% I = zeros(L); 
% I(L+1:L+1:end) = 1;
% I(L,1) = 1;
% vect_order = I*(1:L)';
% a_bar = permute(a_bar,vect_order);
% c = permute(c,vect_order);
%
%--------------------------------------------------------------------------

% Padding
%--------------------------------------------------------------------------
% Mc = size(c);
% Ma = size(a_bar);
% M = max([Mc;Ma]);
% a_temp = a_bar;
% c_temp = c;
% a_bar = zeros(M);
% c = zeros(M);
% for i = 1:length(M)
%     otherdims_a{i} = 1:Ma(i);
%     otherdims_c{i} = 1:Mc(i);
% end
% a_bar(otherdims_a{:}) = a_temp;
% c(otherdims_c{:}) = c_temp;

% Code
%--------------------------------------------------------------------------
Nop_coeff = [0,0,0,-1]; % Cubic S-H 
N = size(a_bar)-1;
mu_k = intval(mu_k_SH(N,lambda)); % \mu_k = lambda - (1-k^2)^2 
[r, Y0,Z0,Z1] = Bounds_int(mu_k,intval(Nop_coeff),intval(a_bar),intval(c),intval(h),intval(e_j),intval(nu));
%--------------------------------------------------------------------------
end

