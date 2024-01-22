function norma = wnorm(a,nu) % return weighted ell-one norm
% WeightedNORM
% a: d-dimensional tensor (d: spatial dimension)
% nu: exponential weights (nu >= 1)
%
if ismatrix(a)
    if isvector(a) % vector case
        n = length(a);
        alp = 2*ones(n,1); alp(1) = 1;
        w = nu.^((0:n-1)').*alp;
    else  % matrix case
        [n1,n2] = size(a); % length of input
        alp = 4*ones(n1,n2); alp(1,:) = 2; alp(:,1) = 2; alp(1,1) = 1;
        v1 = nu.^(0:n1-1); v2 = nu.^(0:n2-1);
        w = (v1'*v2) .* alp;
    end
else % 3D tensor case (TODO: more general code?)
    if exist('intval','file') && isintval(a(1))
        n = size(a); n1 = n(1); n2 = n(2); n3 = n(3);
    else
        [n1,n2,n3] = size(a); % length of input
    end
    alp = 8*ones(n1,n2,n3);
    alp(1,:,:) = 4; alp(:,1,:) = 4; alp(:,:,1) = 4;
    alp(1,1,:) = 2; alp(1,:,1) = 2; alp(:,1,1) = 2;
    alp(1,1,1) = 1;
    v1 = nu.^(0:n1-1); v2 = nu.^(0:n2-1); v3 = nu.^(0:n3-1);
    w = ((permute(v3,[3,1,2]).*ones(n1,n2,n3)) .* (v1'*v2)) .* alp;
end
%
norma = sum(abs(a).*w,"all");