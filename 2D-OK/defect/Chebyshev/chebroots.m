function r = chebroots(c,I)% return roots of Chebyshev series
% c: Two-sided Chebyshev coefficients
% I: domain of the Chebyshev
% 
% NOTE :: available for interval inputs 
% 
arguments
    c
    I = [-1,1];
end
% 
htol = 100*eps;
if c(end)==0
    c(end) = c(end-1);
end
c = -0.5 * c(1:end-1) / c(end);
c(end-1) = c(end-1) + 0.5;

onev = 0.5 * ones(length(c)-1, 1);
A = diag(onev, 1) + diag(onev, -1); % colleague matrix
A(end-1, end) = 1;
if exist('intval','file') && isintval(c(1))
    A = intval(A);
    A(:, 1) = flipud( c );

    [V,D]  = eig(mid(A));
    r = diag(D); % approximated roots

    mask = abs(imag(r)) < htol; % choose reals
    r = real( r(mask) );
    V = V(:,mask);
    mask = abs(r) <= 1 + htol; % Keep roots inside [-1 1]:
    V = V(:,mask);
    % sort the roots
    [r, ind] = sort( r(mask) );
    V = V(:,ind);
    if ( ~isempty(r) )
        r(1) = max(r(1), -1);
        r(end) = min(r(end), 1);
    end
    % Verify roots by verifyeig (require INTLAB)
    n = length(r);
    ir = intval(zeros(n,1));
    for i = 1:n
        ir(i) = verifyeig(A,r(i),V(:,i));
    end
    % Convert roots in [a b]
    a = I(1); b = I(2);
    if ~((a==-1) && (b==1))
        r = (1.0 - ir)*a/2 + (1.0 + ir)*b/2;
    else
        r = ir;
    end
else % non-rigorous numerics
    A(:, 1) = flipud( c );

    r = eig(A);

    mask = abs(imag(r)) < htol;
    r = real( r(mask) );
    % Keep roots inside [-1 1]:
    r = sort( r(abs(r) <= 1 + htol) );
    if ( ~isempty(r) )
        r(1) = max(r(1), -1);
        r(end) = min(r(end), 1);
    end
    % Convert roots in [a b]
    a = I(1); b = I(2);
    if ~((a==-1) && (b==1))
        r = (1.0 - r)*a/2 + (1.0 + r)*b/2;
    end
end
end
