function psi = psi_k(ba_sup,k1,k2,k3,m,nu) % ki starts from index 0
N = size(ba_sup);
N1 = N(1); N2 = N(2); N3 = N(3); 
n1 = N1+k1; n2 = N2+k2; n3 = N3+k3;
%% Index ell
[ell1,ell2,ell3] = ndgrid(1:n1,1:n2,1:n3); % indices F_{k+N} from one (because of picking off nonzero indices)
m1 = m(1); m2 = m(2); m3 = m(3);
ell1(1:m1,1:m2,1:m3) = 0;
ell2(1:m1,1:m2,1:m3) = 0;
ell3(1:m1,1:m2,1:m3) = 0;
ell1 = nonzeros(ell1)-1;
ell2 = nonzeros(ell2)-1;
ell3 = nonzeros(ell3)-1; % F_{k+N}\F_m from zero up to n3-1
% 
%% Index k+ell, k-ell
k1_minus_ell1 = abs(k1-ell1); k1_minus_ell1(k1_minus_ell1>N1) = N1; % up to N1-1. If index is N1, the elemtent must be zero
k1_plus_ell1 = abs(k1+ell1);  k1_plus_ell1(k1_plus_ell1 > N1) = N1;
k2_minus_ell2 = abs(k2-ell2); k2_minus_ell2(k2_minus_ell2>N2) = N2;
k2_plus_ell2 = abs(k2+ell2);  k2_plus_ell2(k2_plus_ell2 > N2) = N2;
k3_minus_ell3 = abs(k3-ell3); k3_minus_ell3(k3_minus_ell3>N3) = N3;
k3_plus_ell3 = abs(k3+ell3);  k3_plus_ell3(k3_plus_ell3 > N3) = N3;
% 
%% weight
% alp = 4*ones(n1,n2); alp(1,:) = 2; alp(:,1) = 2; alp(1,1) = 1;
% v1 = nu.^(0:n1-1); v2 = nu.^(0:n2-1);
% w = (v1'*v2) .* alp;
% w_vec = w(sub2ind([n1,n2],ell1+1,ell2+1));
alp = 8*ones(n1,n2,n3);
alp(1,:,:) = 4; alp(:,1,:) = 4; alp(:,:,1) = 4;
alp(1,1,:) = 2; alp(1,:,1) = 2; alp(:,1,1) = 2;
alp(1,1,1) = 1;
v1 = nu.^(0:n1-1); v2 = nu.^(0:n2-1); v3 = nu.^(0:n3-1);
w = ((permute(v3,[3,1,2]).*ones(n1,n2,n3)) .* (v1'*v2)) .* alp;
w_vec = w(sub2ind([n1,n2,n3],ell1+1,ell2+1,ell3+1));
% 
% %========= computing with chebyshev
% tic
% ba2_ext = intval(zeros(n+1,N1+1,N2+1));
% ba2_ext(1:n,1:N1,1:N2) = ba2;
% % 
% nind = kron((1:n)',ones(length(ell1),1));
% N1ind_minus = kron(ones(n,1),k1_minus_ell1+1);
% N2ind_minus = kron(ones(n,1),k2_minus_ell2+1);
% ba2_k_minus_ell = reshape(ba2_ext(sub2ind(size(ba2_ext),nind,N1ind_minus,N2ind_minus)),length(ell1),n);
% % 
% N1ind_plus = kron(ones(n,1),k1_plus_ell1+1);
% N2ind_plus = kron(ones(n,1),k2_plus_ell2+1);
% ba2_k_plus_ell = reshape(ba2_ext(sub2ind(size(ba2_ext),nind,N1ind_plus,N2ind_plus)),length(ell1),n);
% %  
% psi = max(mag(chebmag((ba2_k_minus_ell + ba2_k_plus_ell)')./w_vec))
% toc
% %===================
% 
% tic
% ba_sup = reshape(chebmag(reshape(ba,n,N1*N2)),N1,N2);
ba_sup_ext = intval(zeros(N1+1,N2+1,N3+1));
ba_sup_ext(1:N1,1:N2,1:N3) = ba_sup;
% ba_sup_k_minus_ell = ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_minus_ell1+1,k2_minus_ell2+1,k3_minus_ell3+1));
% ba_sup_k_plus_ell  = ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_plus_ell1+1,k2_plus_ell2+1,k3_plus_ell3+1));
% psi = max(mag((ba_sup_k_minus_ell+ba_sup_k_plus_ell)./w_vec));
tmp = max([...
  ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_minus_ell1+1,k2_minus_ell2+1,k3_minus_ell3+1)),...
  ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_minus_ell1+1,k2_minus_ell2+1,k3_plus_ell3+1)),...
  ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_minus_ell1+1,k2_plus_ell2+1,k3_minus_ell3+1)),...
  ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_minus_ell1+1,k2_plus_ell2+1,k3_plus_ell3+1)),...
  ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_plus_ell1+1,k2_minus_ell2+1,k3_minus_ell3+1)),...
  ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_plus_ell1+1,k2_minus_ell2+1,k3_plus_ell3+1)),...
  ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_plus_ell1+1,k2_plus_ell2+1,k3_minus_ell3+1)),...
  ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_plus_ell1+1,k2_plus_ell2+1,k3_plus_ell3+1))...
  ], [], 2);
psi = max(mag(tmp./w_vec));
% toc