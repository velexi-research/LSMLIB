%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File:        reinitialization_example_2d.m
% Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:    $Revision$
% Modified:    $Date$
% Description: MATLAB demo program for reinitialization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script demos reinitialization of level set functions to be distance
% functions near the zero level set.
%
% The initial level set function is
%
%    x.^2 + 0.2*y.^2 = 0.04
%
% The boundary condition is a large positive number in both coordinate 
% directions.
%
% In this code, reinitialization is done using WENO5 and TVD RK3. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup environment
clear
format long

% set reinitialization parameters
spatial_derivative_order = 1;
tvdrk_order = 1;
max_iterations = 50;

% set up spatial grid parameters
Nx = 199;
Ny = 199;
ghostcell_width = 3;
Nx_with_ghostcells = Nx+2*ghostcell_width;
Ny_with_ghostcells = Ny+2*ghostcell_width;
x_lo = -1;
x_hi = 1;
y_lo = -1;
y_hi = 1;
dx = (x_hi-x_lo)/Nx;
dy = (y_hi-y_lo)/Ny;
dX = [dx dy];
X = (x_lo-(ghostcell_width-0.5)*dx:dx:x_hi+ghostcell_width*dx)';
Y = (y_lo-(ghostcell_width-0.5)*dy:dy:y_hi+ghostcell_width*dy)';
[y,x] = meshgrid(Y,X);  % reverse order so that x-coordinate has fastest
                        % stride in data arrays


% initialize phi
phi_init = x.^2 + 0.2*y.^2 - 0.04;
phi = phi_init;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reinitialize level set function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fill boundary cells
phi_init(:,1:ghostcell_width) = 1;
phi_init(:,Ny+ghostcell_width+1:end) = 1;
phi_init(1:ghostcell_width,:) = 1;
phi_init(Nx+ghostcell_width+1:end,:) = 1;

phi = reinitializeLevelSetFunction(phi_init, ...
                                   ghostcell_width, ...
                                   dX, ...
                                   max_iterations, ...
                                   spatial_derivative_order, ...
                                   tvdrk_order);

% plot initial level set function 
figure(1); clf;
pcolor(x,y,phi_init);
shading interp
hold on 
contourf(x,y,phi_init,-0.5:0.05:0.5);
contour(x,y,phi_init,[0 0],'m','linewidth',2);
xlabel('x');
ylabel('y');
axis([-1 1 -1 1]);
axis square
color_axis = caxis;  % save color axis for plot of reinitialized phi

% plot reinitialized level set function 
figure(2); clf;
pcolor(x,y,phi);
shading interp
hold on 
contourf(x,y,phi,-0.5:0.05:0.5);
contour(x,y,phi_init,[0 0],'m','linewidth',2);
contour(x,y,phi,[0 0],'c','linewidth',2);
xlabel('x');
xlabel('x');
ylabel('y');
axis([-1 1 -1 1]);
axis square
caxis(color_axis);  % use color axis from plot of initial phi

% plot slices to check if reinitialized level set function is a distance
% function
figure(3); clf;
plot(X, phi(:,ceil(Ny_with_ghostcells/2)));
hold on;
plot(X, phi_init(:,ceil(Ny_with_ghostcells/2)),'r');
xlabel('x');
ylabel('phi');
axis([-1 1 -0.4 0.5]);
title('phi vs. x along line y = 0');

figure(4); clf;
plot(Y,phi(ceil(Nx_with_ghostcells/2),:));
hold on;
plot(Y,phi_init(ceil(Nx_with_ghostcells/2),:),'r');
xlabel('y');
ylabel('phi');
axis([-1 1 -0.4 0.5]);
title('phi vs. y along line x = 0');
