%
% File:        test_computeDistanceFunction2d.m
% Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:    $Revision: 1.7 $
% Modified:    $Date: 2006/08/13 15:45:27 $
% Description: MATLAB test code for computeDistanceFunction2d MEX file
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script tests the computeDistanceFunction2d MATLAB MEX-function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup environment
clear
format long

% plotting parameters
num_plots_per_test = 5;

% grid parameters
Nx = 200;
Ny = 200;
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
% Test Problem 1 
% --------------
%  * interface: two circles centered at (0.25,0.25) and (-0.25,-0.25)
%               both with radius 0.2
%  * velocity:  uniform velocity.  U = (1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Test Problem 1...');
test_problem_num = 1;
[X,Y] = meshgrid(x,y);  
center1 = [0.25 0.25]; radius1 = 0.2;
center2 = [-0.25 -0.25]; radius2 = 0.2;
dist1 = sqrt( (X-center1(1)).^2 + (Y-center1(2)).^2 ) - radius1;
dist2 = sqrt( (X-center2(1)).^2 + (Y-center2(2)).^2 ) - radius2;
phi = min(dist1, dist2);
vel_x = ones(size(X));
vel_y = ones(size(X));
distance_function = computeDistanceFunction2d(phi, dX);

% compare FMM distance function with exact distance function
err = max(max(abs(distance_function-phi)))

% plot results
figure((test_problem_num-1)*num_plots_per_test+1); clf;
contourf(X,Y,distance_function,10);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

figure((test_problem_num-1)*num_plots_per_test+2); clf;
contourf(X,Y,phi,10);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Test Problem 2
% --------------
%  * interface: two circles centered at (0.25,0.25) and (-0.25,-0.25)
%               both with radius 0.2
%  * velocity:  "expanding velocity".  U = (x,y)/r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Test Problem 2...');
test_problem_num = 2;
[X,Y] = meshgrid(x,y);  
center1 = [0.25 0.25]; radius1 = 0.2;
center2 = [-0.25 -0.25]; radius2 = 0.2;
dist1 = sqrt( (X-center1(1)).^2 + (Y-center1(2)).^2 ) - radius1;
dist2 = sqrt( (X-center2(1)).^2 + (Y-center2(2)).^2 ) - radius2;
phi = min(dist1, dist2);
r = sqrt(X.^2 + Y.^2);
idx = find(r ~= 0);
vel_x = zeros(size(X)); vel_y = zeros(size(Y));
vel_x(idx) = X(idx)./r(idx);
vel_y(idx) = Y(idx)./r(idx);
idx = find(r == 0);
vel_x(idx) = 0;
vel_y(idx) = 0;
distance_function = computeDistanceFunction2d(phi, dX);

% compare FMM distance function with exact distance function
err = max(max(abs(distance_function-phi)))

% plot results
figure((test_problem_num-1)*num_plots_per_test+1); clf;
contourf(X,Y,distance_function,20);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

figure((test_problem_num-1)*num_plots_per_test+2); clf;
contourf(X,Y,phi,20);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Test Problem 3
% --------------
%  * interface: two circles centered at (0.25,0.1) and (-0.35,0.25)
%               both with radii of 0.2 and 0.3 respectively
%  * velocity:  "expanding velocity".  U = (x,y)/r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Test Problem 3...');
test_problem_num = 3;
[X,Y] = meshgrid(x,y);  
center1 = [0.25 0.1]; radius1 = 0.2;
center2 = [-0.35 0.25]; radius2 = 0.3;
dist1 = sqrt( (X-center1(1)).^2 + (Y-center1(2)).^2 ) - radius1;
dist2 = sqrt( (X-center2(1)).^2 + (Y-center2(2)).^2 ) - radius2;
phi = min(dist1, dist2);
r = sqrt(X.^2 + Y.^2);
idx = find(r ~= 0);
vel_x = zeros(size(X)); vel_y = zeros(size(Y));
vel_x(idx) = X(idx)./r(idx);
vel_y(idx) = Y(idx)./r(idx);
idx = find(r == 0);
vel_x(idx) = 0;
vel_y(idx) = 0;
distance_function = computeDistanceFunction2d(phi, dX);

% compare FMM distance function with exact distance function
err = max(max(abs(distance_function-phi)))

% plot results
figure((test_problem_num-1)*num_plots_per_test+1); clf;
contourf(X,Y,distance_function,20);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

figure((test_problem_num-1)*num_plots_per_test+2); clf;
contourf(X,Y,phi,20);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

