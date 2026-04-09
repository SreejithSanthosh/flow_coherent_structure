% Example code demonstrating the pipeline to compute the Lagrangian
% deformation for 3D Euclidean flows 

clc; clear; close all

dxyz = 0.1; 
x = 0:dxyz:2*pi; y = 0:dxyz:2*pi; z = 0:dxyz:2*pi;
[X,Y,Z] = ndgrid(x,y,z); [N1,N2,N3] = size(X);

% % Advect tracer particles 
r0 = [X(:);Y(:);Z(:)];
options = odeset('AbsTol', 1e-8, 'RelTol', 1e-6);
[t,r_t] = ode45(@(t,y)velode45(t,y),0:0.1:10,r0,options);
num_steps = length(t); num_ICs = size(r0,1)/2;
r_t = reshape(r_t.',N1,N2,N3,3, num_steps);

%% Compute Lagrangian deformation  

x0 = squeeze(r_t(:,:,:,1,1)); x0 = reshape(x0,[N1,N2,N3]);
y0 = squeeze(r_t(:,:,:,2,1)); y0 = reshape(y0,[N1,N2,N3]);
z0 = squeeze(r_t(:,:,:,3,1)); z0 = reshape(z0,[N1,N2,N3]);

xf = squeeze(r_t(:,:,:,1,end)); xf = reshape(xf,[N1,N2,N3]);
yf = squeeze(r_t(:,:,:,2,end)); yf = reshape(yf,[N1,N2,N3]);
zf = squeeze(r_t(:,:,:,3,end)); zf = reshape(zf,[N1,N2,N3]);

[sv_list,vec0_list] = compute_deform_euclidean3d(x0,y0,z0,xf,yf,zf);

%% Visualize the result
figure('Position',[1583 536 675 547])
fntSz = 24; cr = 2;

% Reshape the FTLE field 
FTLE = log(sv_list(:,:,:,end))/(max(t)-min(t));

scatter3(x0(:),y0(:),z0(:),24,FTLE(:),'filled'); colorbar
xlim([0,2*pi]); ylim([0,2*pi]);zlim([0,2*pi]); axis equal
ax = gca; ax.FontSize = fntSz; ax.LineWidth = 2;
c.FontSize = fntSz; c.LineWidth = 2;
title(sprintf('$$ \\Lambda_{%.1f}^{%.1f}$$',...
    min(t),max(t)),'Interpreter','latex','FontSize',fntSz)

%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%
function dydt = velode45(t,y)
dydt = 0*y;

Nq = size(y,1)/3;
xq = y(1:Nq); yq = y(Nq+1:2*Nq);
zq = y(2*Nq+1:end);

[u,v,w] = velocity(xq,yq,zq,t);

dydt(1:Nq) = u;
dydt(Nq+1:2*Nq) = v;
dydt(2*Nq+1:end) = w;
end 

function [u,v,w] = velocity(X,Y,Z,t)
% Parameters of the ABC flow 
A = sqrt(3); B= sqrt(2); C = 1;

u = A*sin(Z)+C*cos(Y);
v = B*sin(X)+A*cos(Z);
w = C*sin(Y)+B*cos(X);
end 
