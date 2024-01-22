function plot_SH_isosurfaces(a,L)
minval = 1;
%
[N1,N2,N3] = size(a);
L1 = L(1); L2 = L(2); L3 = L(3);
%
pad = 60;
%
N1_pad = pad;
N2_pad = pad;
N3_pad = pad;
%
Nx = N1 + N1_pad;
Ny = N2 + N2_pad;
Nz = N3 + N3_pad;
%
a_pad = zeros(Nx,Ny,Nz);
a_pad(1:N1,1:N2,1:N3) = a;
%
dx = pi/(Nx-1)/L1;
dy = pi/(Ny-1)/L2;
dz = pi/(Nz-1)/L3;
%
[X,Y,Z] = meshgrid(0:Nx-1,0:Ny-1,0:Nz-1);
%
X = dx * X;
Y = dy * Y;
Z = dz * Z;
%
%
aa = zeros(2*Nx-2,2*Ny-2,2*Nz-2);
%
%
Sx = Nx;
Sy = Ny;
Sz = Nz;
%
%
aa(Sx:Sx+Nx-2,Sy:Sy+Ny-2,Sz:Sz+Nz-2)=a_pad(1:end-1,1:end-1,1:end-1);
%
aa(Sx-Nx+1:Sx,Sy:Sy+Ny-2,Sz:Sz+Nz-2)=flip(a_pad(:,1:end-1,1:end-1),1);
aa(Sx:Sx+Nx-2,Sy-Ny+1:Sy,Sz:Sz+Nz-2)=flip(a_pad(1:end-1,:,1:end-1),2);
aa(Sx:Sx+Nx-2,Sy:Sy+Ny-2,Sz-Nz+1:Sz)=flip(a_pad(1:end-1,1:end-1,:),3);
%
aa(Sx-Nx+1:Sx,Sy-Ny+1:Sy,Sz:Sz+Nz-2)=flip(flip(a_pad(:,:,1:end-1),1),2);
aa(Sx-Nx+1:Sx,Sy:Sy+Ny-2,Sz-Nz+1:Sz)=flip(flip(a_pad(:,1:end-1,:),1),3);
aa(Sx:Sx+Nx-2,Sy-Ny+1:Sy,Sz-Nz+1:Sz)=flip(flip(a_pad(1:end-1,:,:),2),3);
%
aa(Sx-Nx+1:Sx,Sy-Ny+1:Sy,Sz-Nz+1:Sz)=flip(flip(flip(a_pad,1),2),3);
u = (2*Nx-2)*(2*Ny-2)*(2*Nz-2)*real(ifftn(ifftshift(aa)));
%
clf
%
data = permute(u(1:Sx,1:Sy,1:Sz),[2,1,3]);
%
% data = smooth3(data,'box',5);
%
isoval = min(data,[],"all")*minval;
% h =
patch(isosurface(X,Y,Z,data,isoval),...
  'FaceColor','blue',...
  'EdgeColor','none',...
  'AmbientStrength',.2,...
  'SpecularStrength',.7,...
  'DiffuseStrength',.4);
% rmpath("/Users/takitoshi/Dropbox/LaptopShare/matlab_work/toolbox/Intlab_V11/gradient/")
% isonormals(data,h)
% addpath("/Users/takitoshi/Dropbox/LaptopShare/matlab_work/toolbox/Intlab_V11/gradient/")
%
patch(isocaps(X,Y,Z,data,isoval),...
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