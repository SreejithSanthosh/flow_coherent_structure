% This is an example of how to use the code 
clear; clc; close all;

% Loading dataset and functions
load('example_data.mat')
addpath('./functions/')
Nt = numel(mesh_time);

% Parameters - need user input
delta = 0.2; % geodesic distance over which to compute deformation
ct0 = 1; % select time-step at which OECS analysis needs to be performed 

% Compute the deformation
mesh_r0 = cell2mat(mesh_r(ct0,:)); % Mesh nodes at ct0
mesh_F0 = cell2mat(mesh_F(ct0)); % Mesh faces at ct0
mesh_v0 = cell2mat(mesh_v(ct0,:)); % vel field at ct0

% Compute the deformation metric
[s2,s1,e2,e1] = compute_strain(mesh_r0,mesh_F0,mesh_v0,delta);

% % Visualize the result
close all; fontSz = 24; camAmp = 8; cam_view = [-33 1.7];

% Define the custom dark background color
bgColor = [13 13 13]/255;

% Set the figure background color
f1 = figure('Position',[72 767 1609 426], 'Color', bgColor);

% --- Subplot 1 ---
ax1 = subplot(1,2,1);
trisurf(mesh_F0,mesh_r0(:,1),mesh_r0(:,2),mesh_r0(:,3),real(s2),'Edgecolor','none')
idx = datasample(1:size(e2,1),1000,'Replace',false); % Undersampling the eigenvectors to display
hold on; 
quiver3(mesh_r0(idx,1),mesh_r0(idx,2),mesh_r0(idx,3),...
    e2(idx,1),e2(idx,2),e2(idx,3),'k','ShowArrowHead','off'); 
hold off
c1 = colorbar; 
c1.Color = 'w'; % Make colorbar text white
axis equal off;
ax1.FontSize = fontSz; 
ax1.Color = bgColor;          % Axes background color
ax1.XColor = 'w'; ax1.YColor = 'w'; ax1.ZColor = 'w'; % Axes text/ticks color
camva(camAmp); shading interp; 
t1 = title(sprintf('$$ s_2(\\mathbf{x},t), \\mathbf{e_2}(\\mathbf{x},t), t = %.1f $$',mesh_time(ct0)),...
    'FontSize',fontSz,'Interpreter','latex'); 
t1.Color = 'w'; % Make title text white
view(cam_view)

% --- Subplot 2 ---
ax2 = subplot(1,2,2);
trisurf(mesh_F0,mesh_r0(:,1),mesh_r0(:,2),mesh_r0(:,3),real(s1),'Edgecolor','none')
idx = datasample(1:size(e1,1),1000,'Replace',false); % Undersampling the eigenvectors to display
hold on; 
quiver3(mesh_r0(idx,1),mesh_r0(idx,2),mesh_r0(idx,3),...
    e1(idx,1),e1(idx,2),e1(idx,3),'k','ShowArrowHead','off'); 
hold off
c2 = colorbar; 
c2.Color = 'w'; 
axis equal off;
ax2.FontSize = fontSz; 
ax2.Color = bgColor;
ax2.XColor = 'w'; ax2.YColor = 'w'; ax2.ZColor = 'w';
camva(camAmp); shading interp; 
t2 = title(sprintf('$$ s_1(\\mathbf{x},t), \\mathbf{e_1}(\\mathbf{x},t), t = %.1f $$',mesh_time(ct0)),...
    'FontSize',fontSz,'Interpreter','latex'); 
t2.Color = 'w';
view(cam_view)