%
% File:        demo_solveEikonalEquation2d.m
% Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:    $Revision: 1.4 $
% Modified:    $Date: 2006/08/13 15:45:27 $
% Description: MATLAB demo code for solveEikonalEquation2d MEX file
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script demos the solveEikonalEquation2d MATLAB MEX-function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup environment
clear
format long

% plotting parameters
num_plots_per_example = 4;

% grid parameters
Nx = 500;
Ny = 500;
x_lo = -1;
x_hi = 1;
y_lo = -1;
y_hi = 1;
dx = (x_hi-x_lo)/Nx;
dy = (y_hi-y_lo)/Ny;
x = (x_lo:dx:x_hi)';
y = (y_lo:dy:y_hi)';
dX = [dx dy];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Example Problem 1 
% -----------------
%  * boundary: circles centered at (0.25,0.25) with radius 0.2
%  * speed:  1 for x < 0 
%            2 for x > 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Example Problem 1...');
demo_problem_num = 1;
[X,Y] = meshgrid(x,y);  
center = [0.25 0.25]; radius = 0.2;
boundary_data = -1*ones(size(X));
idx_bdry = find( abs(sqrt((X-center(1)).^2+(Y-center(2)).^2) - radius) ...
               < max(dX)); 
boundary_data(idx_bdry) = 0; 
boundary_data = single(boundary_data);
speed = single(X);
idx_X_left = find(X<=0);
idx_X_right = find(X>0);
speed(idx_X_left) = 1;
speed(idx_X_right) = 2;
phi = solveEikonalEquation2d(boundary_data, speed, dX);

% plot boundary_data
figure((demo_problem_num-1)*num_plots_per_example+1); clf;
pcolor(X,Y,double(boundary_data));
shading flat;
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

% plot speed function
figure((demo_problem_num-1)*num_plots_per_example+2); clf;
pcolor(X,Y,double(speed));
shading flat;
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

% plot phi
figure((demo_problem_num-1)*num_plots_per_example+3); clf;
pcolor(X,Y,double(phi));
shading interp;
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;


return 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Example Problem 2 
% -----------------
%  * boundary: circle centered at (-1.0,-1.0) with radius 0.2
%  * speed:  0 in box bounded by (0.4,0.1) and (0.6,0.3)
%            0 in box bounded by (-0.6,-0.3) and (-0.4,-0.1)
%            1 otherwise
%  * mask: interior circle centered at (-1.0,-1.0) with radius 0.2
%          region where x-y > 0.75
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Example Problem 2...');
demo_problem_num = 2;
[X,Y] = meshgrid(x,y);  
center = [-1.0 -1.0]; radius = 0.2;
boundary_data = -1*ones(size(X));
idx_bdry = find( abs(sqrt((X-center(1)).^2+(Y-center(2)).^2) - radius) ...
               < max(dX)); 
boundary_data(idx_bdry) = 0; 
boundary_data = single(boundary_data);
speed = single(ones(size(X)));
idx_obs_1 = find( (X<=0.75) & (X>=0.25) & (Y<=0.75) & (Y>=0.25) );
speed(idx_obs_1) = 0.00000;
idx_obs_2 = find( (X>=-0.6) & (X<=-0.3) & (Y>=-0.3) & (Y<=0.0) );
speed(idx_obs_2) = 0.00000;
mask = single(1*ones(size(X)));
idx_mask_1 = find( (sqrt((X-center(1)).^2+(Y-center(2)).^2) - radius) ...
               < -max(dX)); 
mask(idx_mask_1) = -1; 
idx_mask_2 = find( X-Y > 0.75 );
mask(idx_mask_2) = -1; 
phi = solveEikonalEquation2d(boundary_data, speed, dX, mask);

% plot boundary_data
figure((demo_problem_num-1)*num_plots_per_example+1); clf;
pcolor(X,Y,double(boundary_data));
shading flat;
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

% plot speed function
figure((demo_problem_num-1)*num_plots_per_example+2); clf;
pcolor(X,Y,double(speed));
shading flat;
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

% plot mask
figure((demo_problem_num-1)*num_plots_per_example+3); clf;
pcolor(X,Y,double(mask));
shading flat;
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

% plot phi
figure((demo_problem_num-1)*num_plots_per_example+4); clf;
pcolor(X,Y,double(phi));
shading interp;
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
caxis([0 5]);
colorbar;


