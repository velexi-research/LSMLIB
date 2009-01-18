%
% File:        accuracy_solveEikonalEquation2d.m
% Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:    $Revision$
% Modified:    $Date$
% Description: MATLAB test code for solveEikonalEquation2d MEX file
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script checks the order of accuracy for the solveEikonalEquation2d 
% MATLAB MEX-function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup environment
clear
format long

% set zero tolerance
zero_tol = 1000*eps;

% grid sizes
grid_sizes = [50, 100, 200, 400];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example Problem 1
% -----------------
%  * boundary: phi = 0 on circle centered at (0,0) with radius 0.2
%  * speed:    1
%  * mask:     interior of circle
%
% NOTES:
% - Two layer of "boundary values" are provided in order to achieve
%   second-order accuracy in the L-infinity norm for the second-order
%   accurate discretization scheme.
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
  dx = (x_hi-x_lo)/N;
  dy = (y_hi-y_lo)/N;
  x = (x_lo:dx:x_hi)';
  y = (y_lo:dy:y_hi)';
  dX = [dx dy];
  [X,Y] = meshgrid(x,y);  

  % display some information for user ...
  disp_str = sprintf('   solving Example Problem 1 with N = %d', N);
  disp(disp_str);

  % compute exact solution
  phi_exact = sqrt(X.^2 + Y.^2) - radius;
  idx_interior = find( sqrt(X.^2 + Y.^2) - radius < zero_tol);
  phi_exact(idx_interior) = 0;

  boundary_data = -1*ones(size(X));
  idx_bdry = find( abs(sqrt(X.^2 + Y.^2) - radius) < (1+zero_tol)*max(dX)); 
  boundary_data(idx_bdry) = sqrt( X(idx_bdry).^2 + Y(idx_bdry).^2 ) - radius;
  speed = ones(size(boundary_data));
  mask = ones(size(boundary_data));
  mask(idx_interior) = -1;
 
  % solve for phi using first-order discretization
  spatial_discretization_order = 1;
  phi = solveEikonalEquation2d(boundary_data, speed, dX, mask, ...
                               spatial_discretization_order);
  phi(idx_interior) = 0;

  errs1_first_order(i) = max(max(abs(phi-phi_exact)));

  % solve for phi using second-order discretization
  spatial_discretization_order = 2;
  phi = solveEikonalEquation2d(boundary_data, speed, dX, mask, ...
                               spatial_discretization_order);
  phi(idx_interior) = 0;

  errs1_second_order(i) = max(max(abs(phi-phi_exact)));

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

% plot solution at finest grid resolution and points with largest errors
figure(3); clf;
pcolor(X,Y,phi);
shading interp;
axis([0.0 0.25 0.0 0.25]);
pbaspect([1 1 1]);
colorbar;
hold on;
err_second_order_idx = find(errs1_second_order(i) == abs(phi-phi_exact));
plot(X(err_second_order_idx),Y(err_second_order_idx),'go');
plot(X(idx_bdry),Y(idx_bdry),'c+');

% plot error at finest grid resolution and points with largest errors
figure(4); clf;
pcolor(X,Y,phi-phi_exact);
shading interp;
axis([x_lo x_hi y_lo y_hi]);
axis([0.0 0.25 0.0 0.25]);
pbaspect([1 1 1]);
colorbar;
hold on;
err_idx = find(errs1_second_order(i) == abs(phi-phi_exact));
plot(X(err_idx),Y(err_idx),'go');
plot(X(idx_bdry),Y(idx_bdry),'c+');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example Problem 2
% -----------------
%  * boundary: phi = 1 on circle centered at (0,0) with radius 0.7
%  * speed:    1
%  * mask:     exterior of circle centered at (0,0) with radius 0.7
%          and interior of circle centered at (0,0) with radius 0.1
%
% NOTES:
% - Two layer of "boundary values" are provided in order to achieve
%   second-order accuracy in the L-infinity norm for the second-order
%   accurate discretization scheme.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% problem parameters
inner_radius = 0.1;
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
  dx = (x_hi-x_lo)/N;
  dy = (y_hi-y_lo)/N;
  x = (x_lo:dx:x_hi)';
  y = (y_lo:dy:y_hi)';
  dX = [dx dy];
  [X,Y] = meshgrid(x,y);  

  % display some information for user ...
  disp_str = sprintf('   solving Example Problem 2 with N = %d', N);
  disp(disp_str);

  % compute exact solution
  phi_exact = 1 - (sqrt(X.^2 + Y.^2) - outer_radius);
  idx_exterior = find( ...
      (sqrt(X.^2 + Y.^2) - outer_radius > -zero_tol) ...
    | (sqrt(X.^2 + Y.^2) - inner_radius < zero_tol) );
  phi_exact(idx_exterior) = 0;

  boundary_data = -1*ones(size(X));
  idx_bdry = find( ...
      abs(sqrt(X.^2 + Y.^2) - outer_radius) < (1+zero_tol)*max(dX) ); 
  boundary_data(idx_bdry) = 1 ...
    - ( sqrt( X(idx_bdry).^2 + Y(idx_bdry).^2 ) - outer_radius );
  speed = ones(size(boundary_data));
  mask = ones(size(boundary_data));
  mask(idx_exterior) = -1;
 
  % solve for phi using first-order discretization
  spatial_discretization_order = 1;
  phi = solveEikonalEquation2d(boundary_data, speed, dX, mask, ...
                               spatial_discretization_order);
  phi(idx_exterior) = 0;

  errs2_first_order(i) = max(max(abs(phi-phi_exact)));

  % solve for phi using second-order discretization
  spatial_discretization_order = 2;
  phi = solveEikonalEquation2d(boundary_data, speed, dX, mask, ...
                               spatial_discretization_order);
  phi(idx_exterior) = 0;

  errs2_second_order(i) = max(max(abs(phi-phi_exact)));

end

% plot results
figure(5); clf;
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

figure(6); clf;
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

% plot solution at finest grid resolution and points with largest errors
figure(7); clf;
pcolor(X,Y,phi);
shading interp;
axis([0.0 0.7 0.0 0.7]);
pbaspect([1 1 1]);
colorbar;
hold on;
err_second_order_idx = find(errs2_second_order(i) == abs(phi-phi_exact));
plot(X(err_second_order_idx),Y(err_second_order_idx),'go');
plot(X(idx_bdry),Y(idx_bdry),'c+');

% plot error at finest grid resolution and points with largest errors
figure(8); clf;
pcolor(X,Y,phi-phi_exact);
shading interp;
axis([x_lo x_hi y_lo y_hi]);
axis([0.0 0.25 0.0 0.25]);
pbaspect([1 1 1]);
colorbar;
hold on;
err_idx = find(errs2_second_order(i) == abs(phi-phi_exact));
plot(X(err_idx),Y(err_idx),'go');
plot(X(idx_bdry),Y(idx_bdry),'c+');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example Problem 3
% -----------------
%  * boundary: phi = 0 on circle centered at (0,0) with radius 0.2
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
  dx = (x_hi-x_lo)/N;
  dy = (y_hi-y_lo)/N;
  x = (x_lo:dx:x_hi)';
  y = (y_lo:dy:y_hi)';
  dX = [dx dy];
  [X,Y] = meshgrid(x,y);  

  % display some information for user ...
  disp_str = sprintf('   solving Example Problem 3 with N = %d', N);
  disp(disp_str);

  % compute exact solution
  phi_exact = sqrt(X.^2 + Y.^2) - radius;
  idx_interior = find( sqrt(X.^2 + Y.^2) - radius < zero_tol);
  phi_exact(idx_interior) = 0;

  boundary_data = -1*ones(size(X));
  idx_bdry = find( abs(sqrt(X.^2+Y.^2) - radius) < 0.7*max(dX)); 
  boundary_data(idx_bdry) = sqrt( X(idx_bdry).^2 + Y(idx_bdry).^2 ) - radius;
  speed = ones(size(boundary_data));
  mask = ones(size(boundary_data));
  mask(idx_interior) = -1;
 
  % solve for phi using second-order discretization
  spatial_discretization_order = 2;
  phi = solveEikonalEquation2d(boundary_data, speed, dX, mask, ...
                               spatial_discretization_order);
  phi(idx_interior) = 0;

  % compute error in L-infinity and L2 norms
  num_grid_pts = prod(size(X));
  errs3_Linf_norm(i) = max(max(abs(phi-phi_exact)));
  errs3_L2_norm(i) = sqrt( ...
    norm(reshape((phi-phi_exact).^2,1,num_grid_pts),2)*dx*dy );

end

% plot results
figure(9); clf;
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

figure(10); clf;
P_L2 = polyfit(log(grid_sizes),log(errs3_L2_norm),1);
order_L2 = -P_L2(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_L2(1)+P_L2(2)),'k');
hold on;
plot(grid_sizes,errs3_L2_norm,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Example Problem 3: Second-Order Scheme, L_2 Convergence Rate');
order_str = sprintf('Order = %1.1f', order_L2);
text(100,1e-2,order_str);
xlabel('N');
ylabel('L_2 Error');

% plot solution at finest grid resolution and points with largest errors
figure(11); clf;
pcolor(X,Y,phi);
shading interp;
axis([x_lo x_hi y_lo y_hi]);
axis([0.0 0.25 0.0 0.25]);
pbaspect([1 1 1]);
colorbar;
hold on;
err_idx = find(errs3_Linf_norm(i) == abs(phi-phi_exact));
plot(X(err_idx),Y(err_idx),'go');
plot(X(idx_bdry),Y(idx_bdry),'c+');

% plot error at finest grid resolution and points with largest errors
figure(12); clf;
pcolor(X,Y,phi-phi_exact);
shading interp;
axis([x_lo x_hi y_lo y_hi]);
axis([0.0 0.25 0.0 0.25]);
pbaspect([1 1 1]);
colorbar;
hold on;
err_idx = find(errs3_Linf_norm(i) == abs(phi-phi_exact));
plot(X(err_idx),Y(err_idx),'go');
plot(X(idx_bdry),Y(idx_bdry),'c+');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example Problem 4
% -----------------
%  * boundary: phi = 0.0 on circle centered at (0.5,0.5) with radius 0.2
%              phi = 0.2 on circle centered at (-0.5,-0.25) with radius 0.3
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
center1 = [0.5, 0.5];
radius1 = 0.2;
value1  = 0.0;
center2 = [-0.5, -0.25];
radius2 = 0.3;
value2  = 0.2;

% allocate space for errors
errs4_Linf_norm = zeros(size(grid_sizes));
errs4_L2_norm = zeros(size(grid_sizes));

% loop over grid sizes
for i = 1:length(grid_sizes)

  % construct grid 
  N = grid_sizes(i);
  x_lo = -1; x_hi = 1;
  y_lo = -1; y_hi = 1;
  dx = (x_hi-x_lo)/N;
  dy = (y_hi-y_lo)/N;
  x = (x_lo:dx:x_hi)';
  y = (y_lo:dy:y_hi)';
  dX = [dx dy];
  [X,Y] = meshgrid(x,y);  

  % display some information for user ...
  disp_str = sprintf('   solving Example Problem 4 with N = %d', N);
  disp(disp_str);

  % compute indices for boundary and interior
  idx_bdry1 = find( ...
    abs(sqrt((X-center1(1)).^2+(Y-center1(2)).^2) - radius1) < 0.7*max(dX)); 
  idx_bdry2 = find( ...
    abs(sqrt((X-center2(1)).^2+(Y-center2(2)).^2) - radius2) < 0.7*max(dX)); 
  idx_bdry = union(idx_bdry1, idx_bdry2);
  idx_interior1 = ...
     find( sqrt((X-center1(1)).^2 + (Y-center1(2)).^2) - radius1 < zero_tol);
  idx_interior2 = ...
     find( sqrt((X-center2(1)).^2 + (Y-center2(2)).^2) - radius2 < zero_tol);
  idx_interior = union(idx_interior1, idx_interior2);

  % compute exact solution
  phi_exact = min( ...
    sqrt((X-center1(1)).^2 + (Y-center1(2)).^2) - radius1 + value1, ...
    sqrt((X-center2(1)).^2 + (Y-center2(2)).^2) - radius2 + value2);
  phi_exact(idx_interior) = 0;

  boundary_data = -1*ones(size(X));
  boundary_data(idx_bdry1) = ...
      sqrt( (X(idx_bdry1)-center1(1)).^2 + (Y(idx_bdry1)-center1(2)).^2 ) ...
    - radius1 + value1;
  boundary_data(idx_bdry2) = ...
      sqrt( (X(idx_bdry2)-center2(1)).^2 + (Y(idx_bdry2)-center2(2)).^2 ) ...
    - radius2 + value2;
  speed = ones(size(boundary_data));
  mask = ones(size(boundary_data));
  mask(idx_interior) = -1;
 
  % solve for phi using second-order discretization
  spatial_discretization_order = 2;
  phi = solveEikonalEquation2d(boundary_data, speed, dX, mask, ...
                               spatial_discretization_order);
  phi(idx_interior) = 0;

  % compute error in L-infinity and L2 norms
  num_grid_pts = prod(size(X));
  errs4_Linf_norm(i) = max(max(abs(phi-phi_exact)));
  errs4_L2_norm(i) = sqrt( ...
    norm(reshape((phi-phi_exact).^2,1,num_grid_pts),2)*dx*dy );

end

% plot results
figure(13); clf;
P_Linf = polyfit(log(grid_sizes),log(errs4_Linf_norm),1);
order_Linf = -P_Linf(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_Linf(1)+P_Linf(2)),'k');
hold on;
plot(grid_sizes,errs4_Linf_norm,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Example Problem 4: Second-Order Scheme, L_\infty Convergence Rate');
order_str = sprintf('Order = %1.1f', order_Linf);
text(100,0.05,order_str);
xlabel('N');
ylabel('L_\infty Error');

figure(14); clf;
P_L2 = polyfit(log(grid_sizes),log(errs4_L2_norm),1);
order_L2 = -P_L2(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_L2(1)+P_L2(2)),'k');
hold on;
plot(grid_sizes,errs4_L2_norm,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Example Problem 4: Second-Order Scheme, L_2 Convergence Rate');
order_str = sprintf('Order = %1.1f', order_L2);
text(100,1e-2,order_str);
xlabel('N');
ylabel('L_2 Error');

% plot solution at finest grid resolution and points with largest errors
figure(15); clf;
pcolor(X,Y,phi);
shading interp;
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;
hold on;
err_idx = find(errs4_Linf_norm(i) == abs(phi-phi_exact));
plot(X(err_idx),Y(err_idx),'go');
plot(X(idx_bdry),Y(idx_bdry),'c+');

% plot error at finest grid resolution and points with largest errors
figure(16); clf;
pcolor(X,Y,phi-phi_exact);
shading interp;
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;
hold on;
err_idx = find(errs4_Linf_norm(i) == abs(phi-phi_exact));
plot(X(err_idx),Y(err_idx),'go');
plot(X(idx_bdry),Y(idx_bdry),'c+');
