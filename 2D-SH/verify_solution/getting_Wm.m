function [Wm,Wm_J,Wm_at_t1] = getting_Wm(Phi,Psi,r_minus,r_minus_backward,h,n,m,nu)
% Back to two-sided chebyshev
Phi(2:end,:,:) = 2*Phi(2:end,:,:);
Psi(2:end,:,:) = 2*Psi(2:end,:,:);
%
m1m2 = prod(m);
size_of_matrix = m1m2;
%
%%%%%%%%%%% Here should be computed one time %%%%%%%%%%%%%%%
% Creat a mesh based on Chebyshev points of higher order
% tic
divide_num = 2^6; % divide number
theta = linspace(intval('pi'),0,divide_num*n+1); % angles of chebyshev points divided by m
t = 0.5*h*(1+cos(theta)); % values of mesh points
mesh = hull(t(1:end-1),t(2:end)); % construct the mesh
mesh_size = length(mesh);
%
% Chebyshev values on the mesh
i = (1:n)';
cheb_val = (cos(i.*acos(-(h-2*mesh)/h)));
%
% Create indices for mesh
aa = 1:mesh_size; e = ones(1,mesh_size);A = [kron(aa,e);kron(e,aa)];
dom_ind = A(:,A(1,:)>=A(2,:)); % dom_ind(1,:) -> t_index, dom_ind(2,:) -> s_index
%
t_cheb_val = cheb_val(:,dom_ind(1,:));
s_cheb_val = cheb_val(:,dom_ind(2,:));
% toc
%
%%%%%%%%%%%%     Wm bounds     %%%%%%%%%%%%%%
%
num_nodes = length(t_cheb_val);
%
% Phi values on each mesh
% CC0 = (reshape(Phi,n,size_of_matrix^2).')*t_cheb_val + intval(0,r_minus,'midrad');
% Phi_value_on_mesh = reshape(CC0.',num_nodes,size_of_matrix,size_of_matrix);
CC0 = (reshape(Phi,n,size_of_matrix^2).')*t_cheb_val;
Phi_value_on_mesh = reshape(CC0.',num_nodes,size_of_matrix,size_of_matrix);
Phi_value_on_mesh = permute(permute(Phi_value_on_mesh,[2,3,1])+infsup(-1,1)*r_minus,[3,1,2]);

%
% Psi values on each mesh
% CC0_bakckward = (reshape(Psi,n,size_of_matrix^2).')*s_cheb_val + intval(0,r_minus_backward,'midrad');
% Psi_value_on_mesh = reshape(CC0_bakckward.',num_nodes,size_of_matrix,size_of_matrix);
CC0_bakckward = (reshape(Psi,n,size_of_matrix^2).')*s_cheb_val;
Psi_value_on_mesh = reshape(CC0_bakckward.',num_nodes,size_of_matrix,size_of_matrix);
Psi_value_on_mesh = permute(permute(Psi_value_on_mesh,[2,3,1])+infsup(-1,1)*r_minus_backward,[3,1,2]);

%
Phi_Psi_value_on_mesh = intval(zeros(num_nodes,size_of_matrix,size_of_matrix));
% tic
for i=1:m1m2
  for j=1:m1m2
    Phi_Psi_value_on_mesh(:,i,j) = sum(reshape(Phi_value_on_mesh(:,i,:),num_nodes,size_of_matrix).*Psi_value_on_mesh(:,:,j),2);
  end
end
% toc
%
% Wm = max(max(sum(abs(Phi_Psi_value_on_mesh),2),[],3)); % This should be modified in the case of different weights
m1 = m(1); m2 = m(2);
alp = 4*ones(m1,m2); alp(1,:) = 2; alp(:,1) = 2; alp(1,1) = 1;
v1 = nu.^(0:m1-1); v2 = nu.^(0:m2-1);
w_k0 = alp .* (v1'*v2);
w0vec = reshape(w_k0,m1m2,1);
Wm = max(max(permute(sum(permute(abs(Phi_Psi_value_on_mesh),[2,3,1]).*w0vec,1),[2,3,1])./w0vec,[],2));
% 
%%%%%%%%%%%%     Wm_at_t1 bounds     %%%%%%%%%%%%%%
% 
% Phi values on each mesh
% iPhi_at_t1 = sum(intval(Phi),1) + r_minus*infsup(-1,1);
iPhi = intval(Phi);
iPhi(1,:,:) = permute(permute(iPhi(1,:,:),[2,3,1]) + r_minus,[3,1,2]);
iPhi_at_t1 = sum(iPhi,1);
Phi_value_on_mesh = repmat(iPhi_at_t1,num_nodes,1,1);
Wm_at_t1 = max(sum(permute(iPhi_at_t1,[2,3,1]).*w0vec,1)./w0vec');
% 
% Psi values on each mesh
% CC0_bakckward = (reshape(Psi,n,size_of_matrix^2).')*s_cheb_val + r_minus_backward*infsup(-1,1);
% Psi_value_on_mesh = reshape(CC0_bakckward.',num_nodes,size_of_matrix,size_of_matrix);
% 
Phi_Psi_value_on_mesh = intval(zeros(num_nodes,size_of_matrix,size_of_matrix));
% tic
for i=1:m1m2
  for j=1:m1m2
     Phi_Psi_value_on_mesh(:,i,j) = sum(reshape(Phi_value_on_mesh(:,i,:),num_nodes,size_of_matrix).*Psi_value_on_mesh(:,:,j),2);
  end
end
% toc
% 
% Wm_at_t1 = max(max(sum(abs(Phi_Psi_value_on_mesh),2),[],3)); % This should be modified in the case of different weights
Wm_J = max(max(permute(sum(permute(abs(Phi_Psi_value_on_mesh),[2,3,1]).*w0vec,1),[2,3,1])./w0vec,[],2));
% 