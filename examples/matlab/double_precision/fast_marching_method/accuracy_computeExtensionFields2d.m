%
% File:        accuracy_computeExtensionFields2d.m
% Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
%                  Regents of the University of Texas.  All rights reserved.
%              (c) 2009 Kevin T. Chu.  All rights reserved.
% Revision:    $Revision$
% Modified:    $Date$
% Description: MATLAB test code for computeExtensionFields2d MEX file
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script checks the order of accuracy for the 
% computeExtensionFields2d MATLAB MEX-function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup environment
clear
format long

% grid sizes
grid_sizes = [100, 200, 400, 800, 1600];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Example Problem 1 
% -----------------
%  * interface:      circle centered at the origin with radius 0.25.
%  * source_fields:  source1(x,y) = x, source2(x,y) = y
%  * mask:           circle centered at the origin with radius 0.1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% problem parameters
radius = 0.4;
mask_radius = 0.2;
mask_value = -0.5;

% allocate space for errors
errs1_dist_first_order_Linf = zeros(size(grid_sizes));
errs1_dist_first_order_L2 = zeros(size(grid_sizes));
errs1_ext1_first_order_Linf = zeros(size(grid_sizes));
errs1_ext1_first_order_L2 = zeros(size(grid_sizes));
errs1_ext2_first_order_Linf = zeros(size(grid_sizes));
errs1_ext2_first_order_L2 = zeros(size(grid_sizes));

errs1_dist_second_order_Linf = zeros(size(grid_sizes));
errs1_dist_second_order_L2 = zeros(size(grid_sizes));
errs1_ext1_second_order_Linf = zeros(size(grid_sizes));
errs1_ext1_second_order_L2 = zeros(size(grid_sizes));
errs1_ext2_second_order_Linf = zeros(size(grid_sizes));
errs1_ext2_second_order_L2 = zeros(size(grid_sizes));

% loop over grid sizes
for i = 1:length(grid_sizes)

  % grid parameters
  N = grid_sizes(i);
  x_lo = -1;
  x_hi = 1;
  y_lo = -1;
  y_hi = 1;
  dx = (x_hi-x_lo)/N;
  dy = (y_hi-y_lo)/N;
  x = (x_lo:dx:x_hi)';
  y = (y_lo:dy:y_hi)';
  dX = [dx dy];
  [X,Y] = meshgrid(x,y);  

  % set mask
  mask = ones(size(X));
  idx_mask = find( sqrt(X.^2 + Y.^2) - mask_radius < 0);
  mask(idx_mask) = -1; 

  % compute initial level set function
  phi = ( X.^2 + Y.^2 ) - radius^2;

  % compute source fields
  source_fields = cell(2,1);
  source_fields{1} = X;
  source_fields{2} = Y;

  % compute exact distance function
  dist_exact = sqrt( X.^2 + Y.^2 ) - radius; 
  dist_exact(idx_mask) = mask_value;

  % compute solution for extension fields
  ext_field1 = radius*X./(sqrt(X.^2 + Y.^2)); 
  ext_field1(idx_mask) = mask_value;
  ext_field2 = radius*Y./(sqrt(X.^2 + Y.^2));
  ext_field2(idx_mask) = mask_value;

  % solve for distance function using first-order discretization
  spatial_discretization_order = 1;
  [distance_function, extension_fields] = ...
    computeExtensionFields2d(phi, source_fields, dX, mask, ...
                             spatial_discretization_order);
  distance_function(idx_mask) = mask_value;
  extension_fields{1}(idx_mask) = mask_value;
  extension_fields{2}(idx_mask) = mask_value;

  % compute error in L-infinity and L2 norms
  num_grid_pts = prod(size(X));
  err_dist = distance_function-dist_exact;
  err_ext_field1 = extension_fields{1}-ext_field1;
  err_ext_field2 = extension_fields{2}-ext_field2;
  errs1_dist_first_order_Linf(i) = ...
    norm(reshape(err_dist,1,num_grid_pts),inf);
  errs1_dist_first_order_L2(i) = sqrt( ...
    norm(reshape(err_dist.^2,1,num_grid_pts),2)*dx*dy );
  errs1_ext1_first_order_Linf(i) = ...
    norm(reshape(err_ext_field1,1,num_grid_pts),inf);
  errs1_ext1_first_order_L2(i) = sqrt( ...
    norm(reshape(err_ext_field1.^2,1,num_grid_pts),2)*dx*dy );
  errs1_ext2_first_order_Linf(i) = ...
    norm(reshape(err_ext_field2,1,num_grid_pts),inf);
  errs1_ext2_first_order_L2(i) = sqrt( ...
    norm(reshape(err_ext_field2.^2,1,num_grid_pts),2)*dx*dy );

  % solve for distance function using second-order discretization
  spatial_discretization_order = 2;
  [distance_function, extension_fields] = ...
    computeExtensionFields2d(phi, source_fields, dX, mask, ...
                             spatial_discretization_order);
  distance_function(idx_mask) = mask_value;
  extension_fields{1}(idx_mask) = mask_value;
  extension_fields{2}(idx_mask) = mask_value;

  % compute error in L-infinity and L2 norms
  num_grid_pts = prod(size(X));
  err_dist = distance_function-dist_exact;
  err_ext_field1 = extension_fields{1}-ext_field1;
  err_ext_field2 = extension_fields{2}-ext_field2;
  errs1_dist_second_order_Linf(i) = ...
    norm(reshape(err_dist,1,num_grid_pts),inf);
  errs1_dist_second_order_L2(i) = sqrt( ...
    norm(reshape(err_dist.^2,1,num_grid_pts),2)*dx*dy );
  errs1_ext1_second_order_Linf(i) = ...
    norm(reshape(err_ext_field1,1,num_grid_pts),inf);
  errs1_ext1_second_order_L2(i) = sqrt( ...
    norm(reshape(err_ext_field1.^2,1,num_grid_pts),2)*dx*dy );
  errs1_ext2_second_order_Linf(i) = ...
    norm(reshape(err_ext_field2,1,num_grid_pts),inf);
  errs1_ext2_second_order_L2(i) = sqrt( ...
    norm(reshape(err_ext_field2.^2,1,num_grid_pts),2)*dx*dy );

end

% plot results for distance function
figure(1); clf;
P_first_order_Linf = polyfit(log(grid_sizes), ...
                             log(errs1_dist_first_order_Linf),1);
order_first_order_Linf = -P_first_order_Linf(1);
P_first_order_L2 = polyfit(log(grid_sizes),log(errs1_dist_first_order_L2),1);
order_first_order_L2 = -P_first_order_L2(1);
N_plot = [100, 10000];
loglog(N_plot,exp(log(N_plot)*P_first_order_Linf(1)+P_first_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_first_order_L2(1)+P_first_order_L2(2)),'k');
plot(grid_sizes,errs1_dist_first_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_dist_first_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: First-Order Scheme - Distance Function');
order_str = sprintf('Order (L_%s) = %1.2f', '\infty', order_first_order_Linf);
text(300,5e-3,order_str);
order_str = sprintf('Order (L_2) = %1.2f', order_first_order_L2);
text(200,1e-5,order_str);
xlabel('N');
ylabel('Error');

figure(2); clf;
P_second_order_Linf = polyfit(log(grid_sizes), ...
                              log(errs1_dist_second_order_Linf),1);
order_second_order_Linf = -P_second_order_Linf(1);
P_second_order_L2 = polyfit(log(grid_sizes),log(errs1_dist_second_order_L2),1);
order_second_order_L2 = -P_second_order_L2(1);
N_plot = [100, 10000];
loglog(N_plot,exp(log(N_plot)*P_second_order_Linf(1)+P_second_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_second_order_L2(1)+P_second_order_L2(2)),'k');
plot(grid_sizes,errs1_dist_second_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_dist_second_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: Second-Order Scheme - Distance Function');
order_str = sprintf('Order (L_%s) = %1.2f', '\infty', order_second_order_Linf);
text(300,3e-3,order_str);
order_str = sprintf('Order (L_2) = %1.2f', order_second_order_L2);
text(200,1e-6,order_str);
xlabel('N');
ylabel('Error');

% plot error at finest grid resolution and points with largest errors
figure(3); clf;
pcolor(X,Y,abs(err_dist));
shading flat;
axis([x_lo x_hi y_lo y_hi]);
colorbar; 
hold on;
contour(X,Y,distance_function,[0 0],'m');
err_idx = find( ...
  errs1_dist_second_order_Linf(i) == abs(distance_function-dist_exact));
plot(X(err_idx),Y(err_idx),'ko','MarkerFaceColor','k','MarkerSize',13);


% plot results for extension field 1
figure(4); clf;
P_first_order_Linf = polyfit(log(grid_sizes), ...
                             log(errs1_ext1_first_order_Linf),1);
order_first_order_Linf = -P_first_order_Linf(1);
P_first_order_L2 = polyfit(log(grid_sizes),log(errs1_ext1_first_order_L2),1);
order_first_order_L2 = -P_first_order_L2(1);
N_plot = [100, 10000];
loglog(N_plot,exp(log(N_plot)*P_first_order_Linf(1)+P_first_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_first_order_L2(1)+P_first_order_L2(2)),'k');
plot(grid_sizes,errs1_ext1_first_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_ext1_first_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: First-Order Scheme - Extension Field 1');
order_str = sprintf('Order (L_%s) = %1.2f', '\infty', order_first_order_Linf);
text(300,1e-2,order_str);
order_str = sprintf('Order (L_2) = %1.2f', order_first_order_L2);
text(200,1e-5,order_str);
xlabel('N');
ylabel('Error');

figure(5); clf;
P_second_order_Linf = polyfit(log(grid_sizes), ...
                              log(errs1_ext1_second_order_Linf),1);
order_second_order_Linf = -P_second_order_Linf(1);
P_second_order_L2 = polyfit(log(grid_sizes),log(errs1_ext1_second_order_L2),1);
order_second_order_L2 = -P_second_order_L2(1);
N_plot = [100, 10000];
loglog(N_plot,exp(log(N_plot)*P_second_order_Linf(1)+P_second_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_second_order_L2(1)+P_second_order_L2(2)),'k');
plot(grid_sizes,errs1_ext1_second_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_ext1_second_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: Second-Order Scheme - Extension Field 1');
order_str = sprintf('Order (L_%s) = %1.2f', '\infty', order_second_order_Linf);
text(300,4e-3,order_str);
order_str = sprintf('Order (L_2) = %1.2f', order_second_order_L2);
text(200,1e-5,order_str);
xlabel('N');
ylabel('Error');

% plot error at finest grid resolution and points with largest errors
figure(6); clf;
pcolor(X,Y,abs(err_ext_field1));
pcolor(X,Y,err_ext_field1);
shading flat;
axis([x_lo x_hi y_lo y_hi]);
colorbar; 
hold on;
err_idx = find( ...
  errs1_ext1_second_order_Linf(i) == abs(extension_fields{1}-ext_field1));
plot(X(err_idx),Y(err_idx),'ko','MarkerFaceColor','k','MarkerSize',13);


% plot results for extension field 2
figure(7); clf;
P_first_order_Linf = polyfit(log(grid_sizes), ...
                             log(errs1_ext2_first_order_Linf),1);
order_first_order_Linf = -P_first_order_Linf(1);
P_first_order_L2 = polyfit(log(grid_sizes),log(errs1_ext2_first_order_L2),1);
order_first_order_L2 = -P_first_order_L2(1);
N_plot = [100, 10000];
loglog(N_plot,exp(log(N_plot)*P_first_order_Linf(1)+P_first_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_first_order_L2(1)+P_first_order_L2(2)),'k');
plot(grid_sizes,errs1_ext2_first_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_ext2_first_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: First-Order Scheme - Extension Field 2');
order_str = sprintf('Order (L_%s) = %1.2f', '\infty', order_first_order_Linf);
text(300,1e-2,order_str);
order_str = sprintf('Order (L_2) = %1.2f', order_first_order_L2);
text(200,1e-5,order_str);
xlabel('N');
ylabel('Error');

figure(8); clf;
P_second_order_Linf = polyfit(log(grid_sizes), ...
                              log(errs1_ext2_second_order_Linf),1);
order_second_order_Linf = -P_second_order_Linf(1);
P_second_order_L2 = polyfit(log(grid_sizes),log(errs1_ext2_second_order_L2),1);
order_second_order_L2 = -P_second_order_L2(1);
N_plot = [100, 10000];
loglog(N_plot,exp(log(N_plot)*P_second_order_Linf(1)+P_second_order_Linf(2)),'k');
hold on;                       
loglog(N_plot,exp(log(N_plot)*P_second_order_L2(1)+P_second_order_L2(2)),'k');
plot(grid_sizes,errs1_ext2_second_order_Linf,'bo','MarkerSize',14, ...
     'MarkerFaceColor','b');
plot(grid_sizes,errs1_ext2_second_order_L2,'ro','MarkerSize',14, ...
     'MarkerFaceColor','r');
title('Example Problem 1: Second-Order Scheme - Extension Field 2');
order_str = sprintf('Order (L_%s) = %1.2f', '\infty', order_second_order_Linf);
text(300,4e-3,order_str);
order_str = sprintf('Order (L_2) = %1.2f', order_second_order_L2);
text(200,1e-5,order_str);
xlabel('N');
ylabel('Error');

% plot error at finest grid resolution and points with largest errors
figure(9); clf;
pcolor(X,Y,abs(err_ext_field2));
pcolor(X,Y,err_ext_field2);
shading flat;
axis([x_lo x_hi y_lo y_hi]);
colorbar; 
hold on;
err_idx = find( ...
  errs1_ext2_second_order_Linf(i) == abs(extension_fields{2}-ext_field2));
plot(X(err_idx),Y(err_idx),'ko','MarkerFaceColor','k','MarkerSize',13);

