function cinf = getting_cinf(ba_sup,m,nu)
N = size(ba_sup);
N1 = N(1); N2 = N(2); N3 = N(3); 
% padding for extra indices
ba_sup_ext = intval(zeros(N1+1,N2+1,N3+1));
ba_sup_ext(1:N1,1:N2,1:N3) = ba_sup;
% 
%% Index ell
m1 = m(1); m2 = m(2); m3 = m(3);
[ell1,ell2,ell3] = ndgrid(1:m1,1:m2,1:m3); % indices F_{m} from one
ell1(1,1) = 0;
ell2(1,1) = 0;
ell3(1,1) = 0;
ell1 = nonzeros(ell1)-1;
ell2 = nonzeros(ell2)-1; 
ell3 = nonzeros(ell3)-1; % F_m\F_1 from zero up to n3-1
% 
%% Index k
[k1,k2,k3] = ndgrid(1:N1,1:N2,1:N3); % indices F_{N} from one
k1(1:m1,1:m2,1:m3) = 0;
k2(1:m1,1:m2,1:m3) = 0;
k3(1:m1,1:m2,1:m3) = 0; % F_N\F_m
k1 = nonzeros(k1)-1;
k2 = nonzeros(k2)-1;
k3 = nonzeros(k3)-1; % F_N\F_m from zero up to N3-1
% 
%% weight
alp = 8*ones(N1,N2,N3);
alp(1,:,:) = 4; alp(:,1,:) = 4; alp(:,:,1) = 4;
alp(1,1,:) = 2; alp(1,:,1) = 2; alp(:,1,1) = 2;
alp(1,1,1) = 1;
v1 = nu.^(0:N1-1); v2 = nu.^(0:N2-1); v3 = nu.^(0:N3-1);
w = ((permute(v3,[3,1,2]).*ones(N1,N2,N3)) .* (v1'*v2)) .* alp;
% 
%% loop for ell (vectorize w.r.t. k)
tpsi_k = zeros(size(k1));
for l1 = ell1'
  k1_minus_ell1 = abs(k1-l1); k1_minus_ell1(k1_minus_ell1>N1) = N1; % up to N1-1. If index is N1, the elemtent must be zero
  k1_plus_ell1 = abs(k1+l1);  k1_plus_ell1(k1_plus_ell1 > N1) = N1;
  for l2 = ell2'
    k2_minus_ell2 = abs(k2-l2); k2_minus_ell2(k2_minus_ell2>N2) = N2;
    k2_plus_ell2 = abs(k2+l2);  k2_plus_ell2(k2_plus_ell2 > N2) = N2;
    for l3 = ell3'
      k3_minus_ell3 = abs(k3-l3); k3_minus_ell3(k3_minus_ell3>N3) = N3;
      k3_plus_ell3 = abs(k3+l3);  k3_plus_ell3(k3_plus_ell3 > N3) = N3;
      % ba_sup_k_minus_ell = ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_minus_ell1+1,k2_minus_ell2+1,k3_minus_ell3+1));
      % ba_sup_k_plus_ell  = ba_sup_ext(sub2ind([N1+1,N2+1,N3+1],k1_plus_ell1+1,k2_plus_ell2+1,k3_plus_ell3+1));
      % tmp = (ba_sup_k_minus_ell+ba_sup_k_plus_ell)/w(l1+1,l2+1,l3+1);
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
      tmp = tmp/w(l1+1,l2+1,l3+1);
      tpsi_k = max(tpsi_k,tmp);
    end
  end
end
% % compare with the case of ell = (0,0)
ba_sup_k  = ba_sup_ext(sub2ind([N1+1,N2+1],k1+1,k2+1));
tpsi_k = max(tpsi_k,ba_sup_k);
% taking the weighted norm
cinf = 3 * sum(tpsi_k .* w(sub2ind([N1,N2,N3],k1+1,k2+1,k3+1)));
