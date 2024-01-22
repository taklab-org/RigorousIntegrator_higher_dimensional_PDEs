function [ba, delta, a_end] = getting_the_solution_timestepping(N,tspan,a0_init,lambda,L,nu,rigorous)
%%
N1 = N(1); N2 = N(2);
h = tspan(2)-tspan(1);
%
a0 = reshape(a0_init,N1*N2,1);
%
% chebfunpref.setDefaults('factory');
% chebfunpref.setDefaults('fixedLength',n);
% opts = odeset('abstol',1e-16,'reltol',2.22045e-14);
% % u_cheb = chebfun.ode45(@(t,y) F_SH_2D(y,lambda,[N1,N2]),tspan,a0,opts);
% % u_cheb = chebfun.ode113(@(t,y) F_SH_2D(y,lambda,[N1,N2]),tspan,a0,opts);
% u_cheb = chebfun.ode15s(@(t,y) F_SH_2D(y,lambda,[N1,N2]),tspan,a0,opts);
% % u_cheb = chebfun.ode15s(@(t,y) F_SH_3D(y,lambda,[N1,N2,N3]),tspan,a0);
% a = (chebcoeffs(u_cheb)); % (Chebcoeff(including factor 2), Cosine (multi-dimension))
% a_end = reshape(u_cheb(end),N1,N2);
% chebfunpref.setDefaults('factory'); % 3.0235e-09
%
% addpath("defect/Chebyshev/")
opts = odeset('abstol',1e-16,'reltol',2.22045e-14);
f = ode15s(@(t,y) F_SH_2D(y,lambda,N,mid(L)),tspan,a0,opts);
% f = ode113(@(t,y) F_SH_2D(y,lambda,[N1,N2]),tspan,a0,opts);
a = pdecheb2(f,tspan); % Two-sided Chebyshev
n = size(a,1); disp(['n = ',num2str(n)])
a_end = reshape(f.y(:,end).',N1,N2);
% 2.9422e-09(113), 2.9405e-09(89), 3.2228e-09(78), 2.9404e-09(45)
%
if rigorous
  a = intval(a);
else
  L = mid(L);
end
% plotcoeffs(u_cheb),return
% mesh(abs(reshape(u_cheb(end),N1,N2)))
% set(gca, 'ZScale', 'log')
% return
%%%%%%% Compute residual
% ta = [a(1,:);a(2:end,:)/2];
%
% time derivative
rescaleFactork = h/2;
du = chebdiff(a)/rescaleFactork;
du = [du(1,:);du(2:end,:)/2];
du = reshape(du,n-1,N1,N2);
%
if exist('intval','file') && isintval(a(1))
  du_ext = intval(zeros(3*n-2,3*N1-2,3*N2-2));
else
  du_ext = zeros(3*n-2,3*N1-2,3*N2-2);
end
du_ext(1:n-1,1:N1,1:N2) = du;
%
%
ta = [a(1,:);a(2:end,:)/2]; % Two-sided --> One-sided Chebyshev
Fa_ext = F_SH_2D_ext(ta,lambda,n,N,L);
%
defect = du_ext - Fa_ext;
defect(2:end,:,:) = 2*defect(2:end,:,:); % back to Two-sided Chebyshev
delta = wnorm(permute(sum(abs(defect),1),[2,3,1]),nu);
% delta3 = wnorm(reshape(chebmag(reshape(defect,3*n-2,(3*N1-2)*(3*N2-2)),tspan),(3*N1-2),(3*N2-2)),1)
%
ba = reshape(ta,n,N1,N2); % return One-sided Chebyshev
if rigorous
  ba = mid(ba);
end
% rmpath("defect/Chebyshev/");
end