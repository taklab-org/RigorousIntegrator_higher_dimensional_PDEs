close all
clear
clc


% Loading data

%--------------------------------------------------------------------------
load ba_2D_SH.mat
% load ba_2D_SH_alt.mat
% [n,N1,N2] = size(ba);
n_pad = 0;
% ba_ext = (zeros(n+n_pad,N1,N2));
% ba_ext(1:n,:,:) = ba;
% ba = ba_ext;
% Parameters
% -------------------------------------------------------------------------
Nop_coeff = [0,0,0,-1];
h = 2^-7;
nu = 1.4;
lambda = 3;
m = [1,3];
m1m2 = prod(m);
n = size(ba,1);
% L = [1 1 1 1];
% L = [1 1.1];
 L = [1 1];
% Variation Problem
% ------------------------------------------------------------------------
[Phi,r1,Y_01,Z_01,Z_11] = rig_2D_SH_variational(ba,m,Nop_coeff,h,lambda,nu,L);
Phi = reshape(Phi,n+n_pad,m1m2,m1m2);
Phi(2:end,:,:) = 2*Phi(2:end,:,:);


% Adjoint Problem
% ------------------------------------------------------------------------
[Psi,r2,Y_02,Z_02,Z_12] = rig_2D_SH_adjoint(ba,m,Nop_coeff,h,lambda,nu,L);
Psi = reshape(Psi,n+n_pad,m1m2,m1m2);
Psi(2:end,:,:) = 2*Psi(2:end,:,:);
% 
permute(sum(Phi,1),[2 3 1])*permute(sum(Psi,1),[2 3 1]) % this should be identity matrix
norm(eye(m(1)*m(2),m(1)*m(2))-permute(sum(Phi,1),[2 3 1])*permute(sum(Psi,1),[2 3 1]),1)



% Some Plots
% ------------------------------------------------------------------------
figure 
t = linspace(-1,1,100);
norm_vect = zeros(size(t));
ind = 1;
a_plot_full = zeros(size(t));
for k = t
a_plot_full = sum(sum(sum(eval_cheb(permute(ba,[2 3 1]),k))));
test1 = sum(eval_cheb(permute(Phi,[2 3 1]),k),3);
test2 = sum(eval_cheb(permute(Psi,[2 3 1]),k),3);
norm_vect(ind) = norm(eye(m(1)*m(2),m(1)*m(2)) - test1*test2);
ind = ind+1;
end

plot(h/2*t+h/2,norm_vect); hold on
xlabel( ' 0 <= s = t <= h ' )
ylabel( ' || I - \Phi(t) \Psi(s) ||' )
title('m = [2,3]')

M = 50;
figure
[AB,X,Y,T] = plot_ba(ba,linspace(-1,1,M),linspace(0,pi,M),linspace(0,pi,M),L);
x = X(:,:,end); y = Y(:,:,end); ab = AB(:,:,end); surf(x,y,ab); shading interp ;  hold on;
figure
a_cut = ba(:,1:2,1:2);
[AB,X,Y,T] = plot_ba(a_cut,linspace(-1,1,M),linspace(0,pi,M),linspace(0,pi,M),L);
x = X(:,:,end); y = Y(:,:,end); ab = AB(:,:,end); surf(x,y,ab); shading interp ;  hold on;



















