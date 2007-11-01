%
% File:        accuracy_solveEikonalEquation2d.m
% Copyright:   (c) 2005-2008 Kevin T. Chu
% Revision:    $Revision: 1.4 $
% Modified:    $Date: 2006/08/13 15:45:27 $
% Description: MATLAB test code for solveEikonalEquation2d MEX file
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script checks the order of accuracy for the 
% solveEikonalEquation2d MATLAB MEX-function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup environment
clear
format long

% set zero tolerance
zero_tol = 1000*eps;

% grid sizes
grid_sizes = [50, 100, 200, 400, 800];

% allocate space for errors
errs1_first_order = zeros(size(grid_sizes));
errs1_second_order = zeros(size(grid_sizes));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test Problem 1
% --------------
%  * boundary: phi = 0 on circle centered at (0,0) with radius 0.2
%  * speed:    1
%  * mask:     interior of circle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% problem parameters
center = [0.0,0.0]; 
radius = 0.2;

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
  disp_str = sprintf('   solving Test Problem 1 with N = %d', N);
  disp(disp_str);

  % compute exact solution
  phi_exact = sqrt(X.^2 + Y.^2) - radius;
  idx_interior = find( sqrt((X-center(1)).^2+(Y-center(2)).^2) - radius ...
                     < zero_tol);
  phi_exact(idx_interior) = 0;

  boundary_data = -1*ones(size(X));
  idx_bdry = find( abs(sqrt((X-center(1)).^2+(Y-center(2)).^2) - radius) ...
                 < (1+zero_tol)*max(dX)); 
  boundary_data(idx_bdry) = sqrt( (X(idx_bdry)-center(1)).^2 ...
                                + (Y(idx_bdry)-center(2)).^2 ) - radius;
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
title('Test Problem 1: First-Order Scheme');
order_str = sprintf('Order = %1.1f', order_first_order);
text(100,0.05,order_str);
xlabel('N');
ylabel('L_\infty Error');

% plot results
figure(2); clf;
P_second_order = polyfit(log(grid_sizes),log(errs1_second_order),1);
order_second_order = -P_second_order(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_second_order(1)+P_second_order(2)),'k');
hold on;
plot(grid_sizes,errs1_second_order,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Test Problem 1: Second-Order Scheme');
order_str = sprintf('Order = %1.1f', order_second_order);
text(100,1e-2,order_str);
xlabel('N');
ylabel('L_\infty Error');

% plot solution at finest grid resolution and points with largest errors
figure(3); clf;
pcolor(X,Y,phi);
shading interp;
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;
hold on;
err_second_order_idx = find(errs1_second_order(i) == abs(phi-phi_exact));
plot(X(err_second_order_idx),Y(err_second_order_idx),'mo');
%plot(X(idx_interior),Y(idx_interior),'rx');
%plot(X(idx_bdry),Y(idx_bdry),'c+');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test Problem 2
% --------------
%  * boundary: phi = 1 on circle centered at (0,0) with radius 0.7
%  * speed:    1
%  * mask:     exterior of circle centered at (0,0) with radius 0.7
%          and interior of circle centered at (0,0) with radius 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% problem parameters
center = [0.0,0.0]; 
inner_radius = 0.1;
outer_radius = 0.7;

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
  disp_str = sprintf('   solving Test Problem 2 with N = %d', N);
  disp(disp_str);

  % compute exact solution
  phi_exact = 1 - (sqrt(X.^2 + Y.^2) - outer_radius);
  idx_exterior = find( ...
      (sqrt((X-center(1)).^2+(Y-center(2)).^2) - outer_radius > -zero_tol) ...
    | (sqrt((X-center(1)).^2+(Y-center(2)).^2) - inner_radius < zero_tol) );
  phi_exact(idx_exterior) = 0;

  boundary_data = -1*ones(size(X));
  idx_bdry = find( ...
      abs(sqrt((X-center(1)).^2+(Y-center(2)).^2) - outer_radius) ...
    < (1+zero_tol)*max(dX)); 
  boundary_data(idx_bdry) = 1 ...
    - ( sqrt( (X(idx_bdry)-center(1)).^2 + (Y(idx_bdry)-center(2)).^2 ) ...
      - outer_radius );
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
figure(4); clf;
P_first_order = polyfit(log(grid_sizes),log(errs2_first_order),1);
order_first_order = -P_first_order(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_first_order(1)+P_first_order(2)),'k');
hold on;
plot(grid_sizes,errs2_first_order,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Test Problem 2: First-Order Scheme');
order_str = sprintf('Order = %1.1f', order_first_order);
text(100,0.05,order_str);
xlabel('N');
ylabel('L_\infty Error');

% plot results
figure(5); clf;
P_second_order = polyfit(log(grid_sizes),log(errs2_second_order),1);
order_second_order = -P_second_order(1);
N_plot = [10, 1000];
loglog(N_plot,exp(log(N_plot)*P_second_order(1)+P_second_order(2)),'k');
hold on;
plot(grid_sizes,errs2_second_order,'bo','MarkerSize',14,'MarkerFaceColor','b');
title('Test Problem 2: Second-Order Scheme');
order_str = sprintf('Order = %1.1f', order_second_order);
text(100,1e-2,order_str);
xlabel('N');
ylabel('L_\infty Error');

% plot solution at finest grid resolution and points with largest errors
figure(6); clf;
pcolor(X,Y,phi);
shading interp;
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;
hold on;
err_second_order_idx = find(errs2_second_order(i) == abs(phi-phi_exact));
plot(X(err_second_order_idx),Y(err_second_order_idx),'mo');
%plot(X(idx_interior),Y(idx_interior),'rx');
%plot(X(idx_bdry),Y(idx_bdry),'c+');

