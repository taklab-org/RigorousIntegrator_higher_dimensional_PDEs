function fvals = chebvals(c,x,I)% return function values of Chebyshev at x
% c: Two-sided Chebyshev coefficients
% x: evaluation points in I
% I: domain of the Chebyshev
% 
% NOTE :: available for interval inputs, e.g., chebvals(intval(a),0,I))
%         available for multi-dimensional inputs (c & x)
% 
arguments
    c, x
    I = [-1,1];
end
% 
if isempty(x)
    fvals = [];
else
    if exist('intval','file') && isintval(c(1))
        x = intval(x);
    end
a = I(1); b = I(2);
xi = 2*(x-a)/(b-a) - 1;
[M,~] = size(c);
k = 0:M-1;
fvals = cos(k.* acos(xi)) * c;
end