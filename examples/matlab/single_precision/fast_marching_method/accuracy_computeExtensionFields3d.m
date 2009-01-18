%
% File:        accuracy_computeExtensionFields3d.m
% Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:    $Revision$
% Modified:    $Date$
% Description: MATLAB test code for computeExtensionFields3d MEX file
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script checks the order of accuracy for the 
% computeExtensionFields3d MATLAB MEX-function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup environment
clear
format long

% grid sizes
grid_sizes = [100, 150, 200, 250];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Example Problem 1 
% -----------------
%  * interface:      circle centered at the origin with radius 0.25.
%  * source_fields:  source1(x,y,z) = x
%  * mask:           circle centered at the origin with radius 0.2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% problem parameters
radius = 0.5;
mask_radius = 0.2;
mask_value = -0.5;

% allocate space for errors
errs1_dist_first_order_Linf = zeros(size(grid_sizes));
errs1_dist_first_order_L2 = zeros(size(grid_sizes));
errs1_ext1_first_order_Linf = zeros(size(grid_sizes));
errs1_ext1_first_order_L2 = zeros(size(grid_sizes));

errs1_dist_second_order_Linf = zeros(size(grid_sizes));
errs1_dist_second_order_L2 = zeros(size(grid_sizes));
errs1_ext1_second_order_Linf = zeros(size(grid_sizes));
errs1_ext1_second_order_L2 = zeros(size(grid_sizes));

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

  % set mask
  mask = ones(size(X));
  idx_mask = find( sqrt(X.^2 + Y.^2 + Z.^2) - mask_radius < 0);
  mask(idx_mask) = -1; 
  mask = single(mask);

  % compute initial level set function
  phi = single( ( X.^2 + Y.^2  + Z.^2 ) - radius^2 );

  % compute source fields
  source_fields = cell(1,1);
  source_fields{1} = single(X);

  % compute exact distance function
  dist_exact = sqrt( X.^2 + Y.^2 + Z.^2 ) - radius; 
  dist_exact(idx_mask) = mask_value;

  % compute solution for extension fields
  ext_field1 = radius*X./(sqrt(X.^2 + Y.^2 + Z.^2)); 
  ext_field1(idx_mask) = mask_value;

  % display some information for user ...
  disp_str = sprintf('   solving Example Problem 1 with N = %d using first-order method', N);
  disp(disp_str);

  % solve for distance function using first-order discretization
  spatial_discretization_order = 1;
  [distance_function, extension_fields] = ...
    computeExtensionFields3d(phi, source_fields, dX, mask, ...
                             spatial_discretization_order);
  distance_function(idx_mask) = mask_value;
  extension_fields{1}(idx_mask) = mask_value;
  extension_fields{2}(idx_mask) = mask_value;

  % compute error in L-infinity and L2 norms
  num_grid_pts = prod(size(X));
  err_dist = distance_function-dist_exact;
  err_ext_field1 = extension_fields{1}-ext_field1;
  errs1_dist_first_order_Linf(i) = ...
    norm(reshape(err_dist,1,num_grid_pts),inf);
  errs1_dist_first_order_L2(i) = sqrt( ...
    norm(reshape(err_dist.^2,1,num_grid_pts),2)*dx*dy*dz );
  errs1_ext1_first_order_Linf(i) = ...
    norm(reshape(err_ext_field1,1,num_grid_pts),inf);
  errs1_ext1_first_order_L2(i) = sqrt( ...
    norm(reshape(err_ext_field1.^2,1,num_grid_pts),2)*dx*dy*dz );

  % display some information for user ...
  disp_str = sprintf('   solving Example Problem 1 with N = %d using second-order method', N);
  disp(disp_str);

  % solve for distance function using second-order discretization
  spatial_discretization_order = 2;
  [distance_function, extension_fields] = ...
    computeExtensionFields3d(phi, source_fields, dX, mask, ...
                             spatial_discretization_order);
  distance_function(idx_mask) = mask_value;
  extension_fields{1}(idx_mask) = mask_value;
  extension_fields{2}(idx_mask) = mask_value;

  % compute error in L-infinity and L2 norms
  num_grid_pts = prod(size(X));
  err_dist = distance_function-dist_exact;
  err_ext_field1 = extension_fields{1}-ext_field1;
  errs1_dist_second_order_Linf(i) = ...
    norm(reshape(err_dist,1,num_grid_pts),inf);
  errs1_dist_second_order_L2(i) = sqrt( ...
    norm(reshape(err_dist.^2,1,num_grid_pts),2)*dx*dy*dz );
  errs1_ext1_second_order_Linf(i) = ...
    norm(reshape(err_ext_field1,1,num_grid_pts),inf);
  errs1_ext1_second_order_L2(i) = sqrt( ...
    norm(reshape(err_ext_field1.^2,1,num_grid_pts),2)*dx*dy*dz );

end

% plot results for distance function
figure(1); clf;
P_first_order_Linf = polyfit(log(grid_sizes), ...
                             log(errs1_dist_first_order_Linf),1);
order_first_order_Linf = -P_first_order_Linf(1);
P_first_order_L2 = polyfit(log(grid_sizes),log(errs1_dist_first_order_L2),1);
order_first_order_L2 = -P_first_order_L2(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_first_order_Linf(1)+P_first_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_first_order_L2(1)+P_first_order_L2(2)),'k');
plot(grid_sizes,errs1_dist_first_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_dist_first_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: First-Order Scheme - Distance Function');
order_str = sprintf('Order (L_%s) = %1.2f', '\infty', order_first_order_Linf);
text(100,0.05,order_str);
order_str = sprintf('Order (L_2) = %1.2f', order_first_order_L2);
text(20,5e-4,order_str);
xlabel('N');
ylabel('Error');

figure(2); clf;
P_second_order_Linf = polyfit(log(grid_sizes), ...
                              log(errs1_dist_second_order_Linf),1);
order_second_order_Linf = -P_second_order_Linf(1);
P_second_order_L2 = polyfit(log(grid_sizes),log(errs1_dist_second_order_L2),1);
order_second_order_L2 = -P_second_order_L2(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_second_order_Linf(1)+P_second_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_second_order_L2(1)+P_second_order_L2(2)),'k');
plot(grid_sizes,errs1_dist_second_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_dist_second_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: Second-Order Scheme - Distance Function');
order_str = sprintf('Order (L_%s) = %1.2f', '\infty', order_second_order_Linf);
text(100,0.01,order_str);
order_str = sprintf('Order (L_2) = %1.2f', order_second_order_L2);
text(20,1e-4,order_str);
xlabel('N');
ylabel('Error');


% plot results for extension field 1
figure(3); clf;
P_first_order_Linf = polyfit(log(grid_sizes), ...
                             log(errs1_ext1_first_order_Linf),1);
order_first_order_Linf = -P_first_order_Linf(1);
P_first_order_L2 = polyfit(log(grid_sizes),log(errs1_ext1_first_order_L2),1);
order_first_order_L2 = -P_first_order_L2(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_first_order_Linf(1)+P_first_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_first_order_L2(1)+P_first_order_L2(2)),'k');
plot(grid_sizes,errs1_ext1_first_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_ext1_first_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: First-Order Scheme - Extension Field 1');
order_str = sprintf('Order (L_%s) = %1.2f', '\infty', order_first_order_Linf);
text(100,0.05,order_str);
order_str = sprintf('Order (L_2) = %1.2f', order_first_order_L2);
text(20,1e-4,order_str);
xlabel('N');
ylabel('Error');

figure(4); clf;
P_second_order_Linf = polyfit(log(grid_sizes), ...
                              log(errs1_ext1_second_order_Linf),1);
order_second_order_Linf = -P_second_order_Linf(1);
P_second_order_L2 = polyfit(log(grid_sizes),log(errs1_ext1_second_order_L2),1);
order_second_order_L2 = -P_second_order_L2(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_second_order_Linf(1)+P_second_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_second_order_L2(1)+P_second_order_L2(2)),'k');
plot(grid_sizes,errs1_ext1_second_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_ext1_second_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: Second-Order Scheme - Extension Field 1');
order_str = sprintf('Order (L_%s) = %1.2f', '\infty', order_second_order_Linf);
text(100,0.05,order_str);
order_str = sprintf('Order (L_2) = %1.2f', order_second_order_L2);
text(20,1e-4,order_str);
xlabel('N');
ylabel('Error');


