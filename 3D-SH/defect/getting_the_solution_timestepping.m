function [ba, delta, a_end] = getting_the_solution_timestepping(N,tspan,a0_init,lambda,L,nu,rigorous)
%%
N1 = N(1); N2 = N(2); N3 = N(3);
h = tspan(2)-tspan(1);
% 
a0 = reshape(a0_init,N1*N2*N3,1);
% 
% addpath("defect/Chebyshev/")
opts = odeset('abstol',1e-16,'reltol',2.22045e-14);
f = ode15s(@(t,y) F_SH_3D(y,lambda,N,L),tspan,a0,opts);
% f = ode113(@(t,y) F_SH_2D(y,lambda,[N1,N2]),tspan,a0,opts);
% n = 13; a = pdechebcoeffs(f,n,tspan);
a = pdecheb2(f,tspan); % Two-sided Chebyshev
n = size(a,1); disp(['n = ',num2str(n)])
a_end = reshape(f.y(:,end).',N1,N2,N3);
% 2.9422e-09(113), 2.9405e-09(89), 3.2228e-09(78), 2.9404e-09(45)
% 
if rigorous
  a = intval(a); % Switch for rigorous numerics
else
  L = mid(L);
end
% 
%%%%%%% Compute residual
% ta = [a(1,:);a(2:end,:)/2];
% 
% time derivative
rescaleFactork = h/2;
du = chebdiff(a)/rescaleFactork; 
du = [du(1,:);du(2:end,:)/2];
du = reshape(du,n-1,N1,N2,N3);
% 
if exist('intval','file') && isintval(a(1))
  du_ext = intval(zeros(3*n-2,3*N1-2,3*N2-2,3*N3-2));
else
  du_ext = zeros(3*n-2,3*N1-2,3*N2-2,3*N3-2);
end
du_ext(1:n-1,1:N1,1:N2,1:N3) = du;
% 
% 
ta = [a(1,:);a(2:end,:)/2]; % Two-sided --> One-sided Chebyshev
Fa_ext = F_SH_3D_ext(ta,lambda,n,N,L);
% 
defect = du_ext - Fa_ext;
defect(2:end,:,:,:) = 2*defect(2:end,:,:,:); % back to Two-sided Chebyshev
delta = wnorm(permute(sum(abs(defect),1),[2,3,4,1]),nu);
% 
ba = reshape(ta,n,N1,N2,N3); % return One-sided Chebyshev
if rigorous
  ba = mid(ba);
end
end




