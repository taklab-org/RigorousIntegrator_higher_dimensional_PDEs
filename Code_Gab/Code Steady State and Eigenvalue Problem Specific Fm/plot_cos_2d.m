function [Z1,Z2] = plot_cos_2d(a,phi,lambda,L)
x = linspace(0,pi/L(1),100);
y = linspace(0,pi/L(2),100);
[X,Y] = meshgrid(x,y);
N = size(a)-1;
Z1 = zeros(size(X));
Z2 = zeros(size(X));
delta = 4*ones(size(a));
delta(:,1) = 2;
delta(1,:) = 2;
delta(1,1) = 1;
a = a.*delta;
phi = phi.*delta;
for k1 = 0:N(1)
    for k2 = 0:N(2)
        Z1 = Z1 + a(k1+1,k2+1).*cos(k1*L(1)*X).*cos(k2*L(2)*Y);
        Z2 = Z2 + phi(k1+1,k2+1).*cos(k1*L(1)*X).*cos(k2*L(2)*Y);
    end
end

figure
surf(X,Y,Z1);
shading interp
view(2)
axis([0 pi/L(1) 0 pi/L(2)])
figure
surf(X,Y,Z2);
title(['\lambda  = ' num2str(lambda)])
shading interp
view(2)
axis([0 pi/L(1) 0 pi/L(2)])
end

