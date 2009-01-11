%
% File:        accuracy_solveEikonalEquation3d.m
% Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:    $Revision: 1.4 $
% Modified:    $Date: 2006/08/13 15:45:27 $
% Description: MATLAB test code for solveEikonalEquation3d MEX file
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script checks the order of accuracy for the solveEikonalEquation3d 
% MATLAB MEX-function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup environment
clear
format long

% set zero tolerance
zero_tol = 1000*eps;

% grid sizes
grid_sizes = [50, 100, 200];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example Problem 1
% -----------------
%  * boundary: phi = 0 on circle centered at (0,0,0) with radius 0.2
%  * speed:    1
%  * mask:     interior of circle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% problem parameters
radius = 0.2;

% allocate space for errors
errs1_first_order = zeros(size(grid_sizes));
errs1_second_order = zeros(size(grid_sizes));

% loop over grid sizes
for i = 1:length(grid_sizes)

  % construct grid 
  N = grid_sizes(i);
  x_lo = -1; x_hi = 1;
  y_lo = -1; y_hi = 1;
  z_lo = -1; z_hi = 1;
  dx = (x_hi-x_lo)/N;
  dy = (y_hi-y_lo)/N;
  dz = (z_hi-z_lo)/N;
  x = (x_lo:dx:x_hi)';
  y = (y_lo:dy:y_hi)';
  z = (z_lo:dz:z_hi)';
  dX = [dx dy dz];
  [X,Y,Z] = meshgrid(x,y,z);  

  % display some information for user ...
  disp_str = sprintf('   solving Example Problem 1 with N = %d', N);
  disp(disp_str);

  % compute exact solution
  phi_exact = sqrt(X.^2 + Y.^2 + Z.^2) - radius;
  idx_interior = find( sqrt(X.^2 + Y.^2 + Z.^2) - radius < zero_tol);
  phi_exact(idx_interior) = 0;

  boundary_data = -1*ones(size(X));
  idx_bdry = find( ...
      abs(sqrt(X.^2 + Y.^2 + Z.^2) - radius) < (1+zero_tol)*max(dX)); 
  boundary_data(idx_bdry) = sqrt( X(idx_bdry).^2 ...
                                + Y(idx_bdry).^2 ...
                                + Z(idx_bdry).^2 ) - radius;
  boundary_data = single(boundary_data);
  speed = single(ones(size(boundary_data)));
  mask = ones(size(boundary_data));
  mask(idx_interior) = -1;
  mask = single(mask);
 
  % solve for phi using first-order discretization
  spatial_discretization_order = 1;
  phi = solveEikonalEquation3d(boundary_data, speed, dX, mask, ...
                               spatial_discretization_order);
  phi(idx_interior) = 0;

  % compute error for first-order scheme
  num_grid_pts = prod(size(X));
  errs1_first_order(i) = norm(reshape(phi-phi_exact,1,num_grid_pts),inf);

  % solve for phi using second-order discretization
  spatial_discretization_order = 2;
  phi = solveEikonalEquation3d(boundary_data, speed, dX, mask, ...
                               spatial_discretization_order);
  phi(idx_interior) = 0;

  % compute error for second-order scheme
  errs1_second_order(i) = norm(reshape(phi-phi_exact,1,num_grid_pts),inf);

end

% plot results
figure(1); clf;
P_first_order = polyfit(log(grid_sizes),log(errs1_first_order),1);
order_first_order = -P_first_order(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_first_order(1)+P_first_order(2)),'k');
hold on;
plot(grid_sizes,errs1_first_order,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Example Problem 1: First-Order Scheme');
order_str = sprintf('Order = %1.1f', order_first_order);
text(100,0.05,order_str);
xlabel('N');
ylabel('L_\infty Error');

figure(2); clf;
P_second_order = polyfit(log(grid_sizes),log(errs1_second_order),1);
order_second_order = -P_second_order(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_second_order(1)+P_second_order(2)),'k');
hold on;
plot(grid_sizes,errs1_second_order,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Example Problem 1: Second-Order Scheme');
order_str = sprintf('Order = %1.1f', order_second_order);
text(100,1e-2,order_str);
xlabel('N');
ylabel('L_\infty Error');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example Problem 2
% -----------------
%  * boundary: phi = 1 on circle centered at (0,0,0) with radius 0.7
%  * speed:    1
%  * mask:     exterior of circle centered at (0,0,0) with radius 0.7
%          and interior of circle centered at (0,0,0) with radius 0.25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% problem parameters
inner_radius = 0.25;
outer_radius = 0.7;

% allocate space for errors
errs2_first_order = zeros(size(grid_sizes));
errs2_second_order = zeros(size(grid_sizes));

% loop over grid sizes
for i = 1:length(grid_sizes)

  % construct grid 
  N = grid_sizes(i);
  x_lo = -1; x_hi = 1;
  y_lo = -1; y_hi = 1;
  z_lo = -1; z_hi = 1;
  dx = (x_hi-x_lo)/N;
  dy = (y_hi-y_lo)/N;
  dz = (z_hi-z_lo)/N;
  x = (x_lo:dx:x_hi)';
  y = (y_lo:dy:y_hi)';
  z = (z_lo:dz:z_hi)';
  dX = [dx dy dz];
  [X,Y,Z] = meshgrid(x,y,z);  

  % display some information for user ...
  disp_str = sprintf('   solving Example Problem 2 with N = %d', N);
  disp(disp_str);

  % compute exact solution
  phi_exact = 1 - (sqrt(X.^2 + Y.^2 + Z.^2) - outer_radius);
  idx_exterior = find( ...
      ( sqrt(X.^2 + Y.^2 + Z.^2) - outer_radius > -zero_tol ) ...
    | ( sqrt(X.^2 + Y.^2 + Z.^2) - inner_radius < zero_tol ) );
  phi_exact(idx_exterior) = 0;

  boundary_data = -1*ones(size(X));
  idx_bdry = find( ...
      abs( sqrt(X.^2 + Y.^2 + Z.^2) - outer_radius ) < (1+zero_tol)*max(dX)); 
  boundary_data(idx_bdry) = 1 - ( sqrt( X(idx_bdry).^2 ...
                                      + Y(idx_bdry).^2 ...
                                      + Z(idx_bdry).^2 ) - outer_radius);
  boundary_data = single(boundary_data);
  speed = single(ones(size(boundary_data)));
  mask = ones(size(boundary_data));
  mask(idx_exterior) = -1;
  mask = single(mask);
 
  % solve for phi using first-order discretization
  spatial_discretization_order = 1;
  phi = solveEikonalEquation3d(boundary_data, speed, dX, mask, ...
                               spatial_discretization_order);
  phi(idx_exterior) = 0;

  % compute error for first-order scheme
  num_grid_pts = prod(size(X));
  errs2_first_order(i) = norm(reshape(phi-phi_exact,1,num_grid_pts),inf);

  % solve for phi using second-order discretization
  spatial_discretization_order = 2;
  phi = solveEikonalEquation3d(boundary_data, speed, dX, mask, ...
                               spatial_discretization_order);
  phi(idx_exterior) = 0;

  % compute error for second-order scheme
  errs2_second_order(i) = norm(reshape(phi-phi_exact,1,num_grid_pts),inf);

end

% plot results
figure(3); clf;
P_first_order = polyfit(log(grid_sizes),log(errs2_first_order),1);
order_first_order = -P_first_order(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_first_order(1)+P_first_order(2)),'k');
hold on;
plot(grid_sizes,errs2_first_order,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Example Problem 2: First-Order Scheme');
order_str = sprintf('Order = %1.1f', order_first_order);
text(100,0.05,order_str);
xlabel('N');
ylabel('L_\infty Error');

figure(4); clf;
P_second_order = polyfit(log(grid_sizes),log(errs2_second_order),1);
order_second_order = -P_second_order(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_second_order(1)+P_second_order(2)),'k');
hold on;
plot(grid_sizes,errs2_second_order,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Example Problem 2: Second-Order Scheme');
order_str = sprintf('Order = %1.1f', order_second_order);
text(100,1e-2,order_str);
xlabel('N');
ylabel('L_\infty Error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example Problem 3
% -----------------
%  * boundary: phi = 0 on circle centered at (0,0,0) with radius 0.2
%  * speed:    1
%  * mask:     interior of circle
%
% NOTES:
% - Only one layer of boundary values is provided at some points on
%   the boundary to demonstrate that the second-order accurate
%   discretization scheme drops to first-order accuracy in the
%   L-infinity norm.  The solution remains second-order accurate
%   in the L2 norm.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% problem parameters
radius = 0.2;

% allocate space for errors
errs3_Linf_norm = zeros(size(grid_sizes));
errs3_L2_norm = zeros(size(grid_sizes));

% loop over grid sizes
for i = 1:length(grid_sizes)

  % construct grid
  N = grid_sizes(i);
  x_lo = -1; x_hi = 1;
  y_lo = -1; y_hi = 1;
  z_lo = -1; z_hi = 1;
  dx = (x_hi-x_lo)/N;
  dy = (y_hi-y_lo)/N;
  dz = (z_hi-z_lo)/N;
  x = (x_lo:dx:x_hi)';
  y = (y_lo:dy:y_hi)';
  z = (z_lo:dz:z_hi)';
  dX = [dx dy dz];
  [X,Y,Z] = meshgrid(x,y,z);

  % display some information for user ...
  disp_str = sprintf('   solving Example Problem 3 with N = %d', N);
  disp(disp_str);

  % compute exact solution
  phi_exact = sqrt(X.^2 + Y.^2 + Z.^2) - radius;
  idx_interior = find( sqrt(X.^2 + Y.^2 + Z.^2) - radius < zero_tol);
  phi_exact(idx_interior) = 0;

  boundary_data = -1*ones(size(X));
  idx_bdry = find( abs(sqrt(X.^2 + Y.^2 + Z.^2) - radius) < 0.75*max(dX));
  boundary_data(idx_bdry) = sqrt( X(idx_bdry).^2 ...
                                + Y(idx_bdry).^2 ...
                                + Z(idx_bdry).^2 ) - radius;
  boundary_data = single(boundary_data);
  speed = single(ones(size(boundary_data)));
  mask = ones(size(boundary_data));
  mask(idx_interior) = -1;
  mask = single(mask);

  % solve for phi using second-order discretization
  spatial_discretization_order = 2;
  phi = solveEikonalEquation3d(boundary_data, speed, dX, mask, ...
                               spatial_discretization_order);
  phi(idx_interior) = 0;

  % compute error in L-infinity and L2 norms
  num_grid_pts = prod(size(X));
  errs3_Linf_norm(i) = norm(reshape(phi-phi_exact,1,num_grid_pts),inf);
  errs3_L2_norm(i) = sqrt( ...
    norm(reshape((phi-phi_exact).^2,1,num_grid_pts),2)*dx*dy*dz );

end

% plot results
figure(5); clf;
P_Linf = polyfit(log(grid_sizes),log(errs3_Linf_norm),1);
order_Linf = -P_Linf(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_Linf(1)+P_Linf(2)),'k');
hold on;
plot(grid_sizes,errs3_Linf_norm,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Example Problem 3: Second-Order Scheme, L_\infty Convergence Rate');
order_str = sprintf('Order = %1.1f', order_Linf);
text(100,0.05,order_str);
xlabel('N');
ylabel('L_\infty Error');

figure(6); clf;
P_L2 = polyfit(log(grid_sizes),log(errs3_L2_norm),1);
order_L2 = -P_L2(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_L2(1)+P_L2(2)),'k');
hold on;
plot(grid_sizes,errs3_L2_norm,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Example Problem 3: Second-Order Scheme, L_2 Convergence Rate');
order_str = sprintf('Order = %1.1f', order_L2);
text(100,1e-3,order_str);
xlabel('N');
ylabel('L_2 Error');

