%
% File:        accuracy_computeDistanceFunction3d.m
% Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:    $Revision: 1.7 $
% Modified:    $Date: 2006/08/13 15:45:27 $
% Description: MATLAB test code for computeDistanceFunction3d MEX file
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script checks the order of accuracy for the 
% computeDistanceFunction3d MATLAB MEX-function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup environment
clear
format long

% grid sizes
grid_sizes = [25, 50, 100, 200];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Example Problem 1 
% -----------------
%  * interface: circle centered at origin with radius 0.4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% problem parameters
center = [0.0 0.0 0.0]; radius = 0.4;

% allocate space for errors
errs1_first_order_Linf = zeros(size(grid_sizes));
errs1_first_order_L2 = zeros(size(grid_sizes));
errs1_second_order_Linf = zeros(size(grid_sizes));
errs1_second_order_L2 = zeros(size(grid_sizes));

% loop over grid sizes
for i = 1:length(grid_sizes)

  % grid parameters
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

  % compute initial level set function
  phi = ( (X-center(1)).^2 + (Y-center(2)).^2 + (Z-center(3)).^2 ) - radius^2;
  phi = single(phi);

  % compute exact distance function
  dist_exact = ...
    sqrt( (X-center(1)).^2 + (Y-center(2)).^2 + (Z-center(3)).^2 ) - radius;

  % solve for distance function using first-order discretization
  mask = [];
  spatial_discretization_order = 1;
  distance_function = computeDistanceFunction3d(phi, dX, mask, ...
                        spatial_discretization_order);

  % compute error in L-infinity and L2 norms
  err = distance_function-dist_exact;
  num_grid_pts = prod(size(X));
  errs1_first_order_Linf(i) = ...
    norm(reshape(err,1,num_grid_pts),inf);
  errs1_first_order_L2(i) = sqrt( ...
    norm(reshape(err.^2,1,num_grid_pts),2)*dx*dy*dz );

  % solve for distance function using second-order discretization
  mask = [];
  spatial_discretization_order = 2;
  distance_function = computeDistanceFunction3d(phi, dX, mask, ...
                        spatial_discretization_order);

  % compute error in L-infinity and L2 norms
  err = distance_function-dist_exact;
  num_grid_pts = prod(size(X));
  errs1_second_order_Linf(i) = ...
    norm(reshape(err,1,num_grid_pts),inf);
  errs1_second_order_L2(i) = sqrt( ...
    norm(reshape(err.^2,1,num_grid_pts),2)*dx*dy*dz );

end

% plot results
figure(1); clf;
P_first_order_Linf = polyfit(log(grid_sizes),log(errs1_first_order_Linf),1);
order_first_order_Linf = -P_first_order_Linf(1);
P_first_order_L2 = polyfit(log(grid_sizes),log(errs1_first_order_L2),1);
order_first_order_L2 = -P_first_order_L2(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_first_order_Linf(1)+P_first_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_first_order_L2(1)+P_first_order_L2(2)),'k');
plot(grid_sizes,errs1_first_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_first_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: First-Order Scheme');
order_str = sprintf('Order (L_%s) = %1.1f', '\infty', order_first_order_Linf);
text(100,0.05,order_str);
order_str = sprintf('Order (L_2) = %1.1f', order_first_order_L2);
text(20,7e-4,order_str);
xlabel('N');
ylabel('Error');

% plot results
figure(2); clf;
P_second_order_Linf = polyfit(log(grid_sizes),log(errs1_second_order_Linf),1);
order_second_order_Linf = -P_second_order_Linf(1);
P_second_order_L2 = polyfit(log(grid_sizes),log(errs1_second_order_L2),1);
order_second_order_L2 = -P_second_order_L2(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_second_order_Linf(1)+P_second_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_second_order_L2(1)+P_second_order_L2(2)),'k');
plot(grid_sizes,errs1_second_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_second_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: Second-Order Scheme');
order_str = sprintf('Order (L_%s) = %1.1f', '\infty', order_second_order_Linf);
text(100,0.01,order_str);
order_str = sprintf('Order (L_2) = %1.1f', order_second_order_L2);
text(20,1e-4,order_str);
xlabel('N');
ylabel('Error');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Example Problem 2
% -----------------
%  * interface: two circles centered at (0.25,0.25,0.25) and 
%               (-0.25,-0.25,-0.25) both with radius 0.2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% problem parameters
center1 = [0.25 0.25 0.25]; radius1 = 0.2;
center2 = [-0.25 -0.25 -0.25]; radius2 = 0.3;

% allocate space for errors
errs2_first_order_Linf = zeros(size(grid_sizes));
errs2_first_order_L2 = zeros(size(grid_sizes));
errs2_second_order_Linf = zeros(size(grid_sizes));
errs2_second_order_L2 = zeros(size(grid_sizes));

% loop over grid sizes
for i = 1:length(grid_sizes)

  % grid parameters
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

  % compute initial level set function
  dist1 = ( (X-center1(1)).^2 + (Y-center1(2)).^2 + (Z-center1(3)).^2 ) ...
        - radius1^2;
  dist2 = ( (X-center2(1)).^2 + (Y-center2(2)).^2 + (Z-center2(3)).^2 ) ...
        - radius2^2;
  phi = single( min(dist1, dist2) );

  % compute exact distance function
  dist1 = sqrt( (X-center1(1)).^2 + (Y-center1(2)).^2 + (Z-center1(3)).^2 ) ...
        - radius1;
  dist2 = sqrt( (X-center2(1)).^2 + (Y-center2(2)).^2 + (Z-center2(3)).^2 ) ...
        - radius2;
  dist_exact = min(dist1, dist2);

  % solve for distance function using first-order discretization
  mask = [];
  spatial_discretization_order = 1;
  distance_function = computeDistanceFunction3d(phi, dX, mask, ...
                        spatial_discretization_order);

  % compute error in L-infinity and L2 norms
  err = distance_function-dist_exact;
  num_grid_pts = prod(size(X));
  errs2_first_order_Linf(i) = ...
    norm(reshape(err,1,num_grid_pts),inf);
  errs2_first_order_L2(i) = sqrt( ...
    norm(reshape(err.^2,1,num_grid_pts),2)*dx*dy*dz );

  % solve for distance function using second-order discretization
  mask = [];
  spatial_discretization_order = 2;
  distance_function = computeDistanceFunction3d(phi, dX, mask, ...
                        spatial_discretization_order);

  % compute error in L-infinity and L2 norms
  err = distance_function-dist_exact;
  num_grid_pts = prod(size(X));
  errs2_second_order_Linf(i) = ...
    norm(reshape(err,1,num_grid_pts),inf);
  errs2_second_order_L2(i) = sqrt( ...
    norm(reshape(err.^2,1,num_grid_pts),2)*dx*dy*dz );

end

% plot results
figure(3); clf;
P_first_order_Linf = polyfit(log(grid_sizes),log(errs2_first_order_Linf),1);
order_first_order_Linf = -P_first_order_Linf(1);
P_first_order_L2 = polyfit(log(grid_sizes),log(errs2_first_order_L2),1);
order_first_order_L2 = -P_first_order_L2(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_first_order_Linf(1)+P_first_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_first_order_L2(1)+P_first_order_L2(2)),'k');
plot(grid_sizes,errs2_first_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs2_first_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 2: First-Order Scheme');
order_str = sprintf('Order (L_%s) = %1.1f', '\infty', order_first_order_Linf);
text(100,0.05,order_str);
order_str = sprintf('Order (L_2) = %1.1f', order_first_order_L2);
text(20,7e-4,order_str);
xlabel('N');
ylabel('Error');

% plot results
figure(4); clf;
P_second_order_Linf = polyfit(log(grid_sizes),log(errs2_second_order_Linf),1);
order_second_order_Linf = -P_second_order_Linf(1);
P_second_order_L2 = polyfit(log(grid_sizes),log(errs2_second_order_L2),1);
order_second_order_L2 = -P_second_order_L2(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_second_order_Linf(1)+P_second_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_second_order_L2(1)+P_second_order_L2(2)),'k');
plot(grid_sizes,errs2_second_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs2_second_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 2: Second-Order Scheme');
order_str = sprintf('Order (L_%s) = %1.1f', '\infty', order_second_order_Linf);
text(100,0.01,order_str);
order_str = sprintf('Order (L_2) = %1.1f', order_second_order_L2);
text(20,1e-4,order_str);
xlabel('N');
ylabel('Error');

