function [] = plot_cos_3d(a,phi,lambda,L)
x = linspace(0,pi/L(1),100);
y = linspace(0,pi/L(2),100);
z = linspace(0,pi/L(3),100);
[Y,X,Z] = meshgrid(x,y,z);
N = size(a)-1;
Z1 = zeros(size(X));
Z2 = zeros(size(X));

delta = 8*ones(size(a));

delta(:,:,1) = 4;
delta(:,1,:) = 4;
delta(1,:,:) = 4;

delta(:,1,1) = 2;
delta(1,1,:) = 2;
delta(1,:,1) = 2;

delta(1,1,1) = 1;

a = a.*delta;
phi = phi.*delta;

for k1 = 0:N(1)
    for k2 = 0:N(2)
        for k3 = 0:N(3)
            Z1 = Z1 + a(k1+1,k2+1,k3+1).*cos(k1*L(1)*X).*cos(k2*L(2)*Y).*cos(k3*L(3)*Z);
            Z2 = Z2 + phi(k1+1,k2+1,k3+1).*cos(k1*L(1)*X).*cos(k2*L(2)*Y).*cos(k3*L(3)*Z);
        end
    end
end
figure 
minval = 1;
isoval = min(Z1,[],"all")*minval;
patch(isosurface(X,Y,Z,Z1,isoval),...
  'FaceColor','blue',...
  'EdgeColor','none',...
  'AmbientStrength',.2,...
  'SpecularStrength',.7,...
  'DiffuseStrength',.4);
% rmpath("/Users/takitoshi/Dropbox/LaptopShare/matlab_work/toolbox/Intlab_V11/gradient/")
% isonormals(data,h)
% addpath("/Users/takitoshi/Dropbox/LaptopShare/matlab_work/toolbox/Intlab_V11/gradient/")
%
patch(isocaps(X,Y,Z,Z1,isoval),...
  'FaceColor','interp',...
  'EdgeColor','none')
colormap jet
% colormap default
% daspect([1,1,1])
axis tight
view(3)
% % % camlight right
% % camlight left
% camlight headlight
%
camlight
camlight(-80,-10)
lighting gouraud
%
xlabel('$x_1$','interpreter','latex', 'FontSize', 30)
ylabel('$x_2$','interpreter', 'latex', 'FontSize', 30)
zlabel('$x_3$','interpreter', 'latex', 'FontSize', 30)
set(gca,'FontSize',20)
colorbar

figure 
minval = 1;
isoval = min(Z2,[],"all")*minval;
patch(isosurface(X,Y,Z,Z2,isoval),...
  'FaceColor','blue',...
  'EdgeColor','none',...
  'AmbientStrength',.2,...
  'SpecularStrength',.7,...
  'DiffuseStrength',.4);
% rmpath("/Users/takitoshi/Dropbox/LaptopShare/matlab_work/toolbox/Intlab_V11/gradient/")
% isonormals(data,h)
% addpath("/Users/takitoshi/Dropbox/LaptopShare/matlab_work/toolbox/Intlab_V11/gradient/")
%
patch(isocaps(X,Y,Z,Z2,isoval),...
  'FaceColor','interp',...
  'EdgeColor','none')
colormap jet
% colormap default
% daspect([1,1,1])
axis tight
view(3)
% % % camlight right
% % camlight left
% camlight headlight
%
camlight
camlight(-80,-10)
lighting gouraud
%
xlabel('$x_1$','interpreter','latex', 'FontSize', 30)
ylabel('$x_2$','interpreter', 'latex', 'FontSize', 30)
zlabel('$x_3$','interpreter', 'latex', 'FontSize', 30)
set(gca,'FontSize',20)
colorbar
title(['\lambda  = ' num2str(lambda)])
end

