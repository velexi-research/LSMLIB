%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File:        const_adv_2d_gaussian.m
% Copyright:   (c) 2005-2006 Kevin T. Chu 
% Revision:    $Revision: 1.4 $
% Modified:    $Date: 2006/04/22 12:02:54 $
% Description: MATLAB test program for 2D ENO/WENO spatial derivatives 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script tests the 2D ENO/WENO derivative calculations using 
% the advection equation:
%
%    u_t + v_x u_x + v_y u_y = 0
%
% where (v_x,v_y) is a constant advection velocity.  The initial condition
% is Gaussian bump:
%
%   u(x,y,t=0) = 1/alpha*exp( -alpha*(x^2+y^2) )
%
% and the boundary conditions are periodic in both coordinate directions.
%
% In this code, time advection is done using forward euler (TVD RK1), 
% TVD RK2, and TVD RK3.
%
% Kevin Chu
% Dept of Mathematics, MIT
% April 2005
%

% setup environment
clear
format long

% set up physical parameters for velocity and width of initial conditions
V_x = 0.33; V_y = 0;
V_x = 0; V_y = 0.33;
V_x = 0.33; V_y = -0.48;
alpha = 50;

% set up spatial grid parameters
Nx = 50;
Ny = 50;
max_ENO_order = 3;
ghostcell_width = max_ENO_order;
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

% set advection velocity function
v_x = V_x*ones(Nx_with_ghostcells,Ny_with_ghostcells);
v_y = V_y*ones(Nx_with_ghostcells,Ny_with_ghostcells);

% set up time integration parameters
cfl_RK1 = 0.4;
cfl_RK2 = 0.6;
cfl_RK3 = 0.6;
dt_RK1 = cfl_RK1/(abs(V_x)/dx+abs(V_y)/dy);
dt_RK2 = cfl_RK2/(abs(V_x)/dx+abs(V_y)/dy);
dt_RK3 = cfl_RK3/(abs(V_x)/dx+abs(V_y)/dy);
t_i = 0;
t_f = 5;
t_RK1 = t_i:dt_RK1:t_f;
if (t_RK1(end) ~= t_f)
  t_RK1 = [t_RK1 t_f];
end
t_RK2 = t_i:dt_RK2:t_f;
if (t_RK2(end) ~= t_f)
  t_RK2 = [t_RK2 t_f];
end
t_RK3 = t_i:dt_RK3:t_f;
if (t_RK3(end) ~= t_f)
  t_RK3 = [t_RK3 t_f];
end


% initialize u
u_init = 1/alpha*exp(-alpha*(x.^2+y.^2));
u_init = 1/alpha*exp(-alpha*(x.^2+y.^2));
u_ENO1_2d_RK1 = u_init;
u_ENO2_2d_RK1 = u_ENO1_2d_RK1;
u_ENO3_2d_RK1 = u_ENO1_2d_RK1;
u_WENO5_2d_RK1 = u_ENO1_2d_RK1;
u_ENO1_2d_RK2 = u_ENO1_2d_RK1;
u_ENO2_2d_RK2 = u_ENO1_2d_RK1;
u_ENO3_2d_RK2 = u_ENO1_2d_RK1;
u_WENO5_2d_RK2 = u_ENO1_2d_RK1;
u_ENO1_2d_RK3 = u_ENO1_2d_RK1;
u_ENO2_2d_RK3 = u_ENO1_2d_RK1;
u_ENO3_2d_RK3 = u_ENO1_2d_RK1;
u_WENO5_2d_RK3 = u_ENO1_2d_RK1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for forward euler (TVD RK1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t_RK1)

  % fill boundary cells
  u_ENO1_2d_RK1(:,1:ghostcell_width) = ...
    u_ENO1_2d_RK1(:,Ny+1:ghostcell_width+Ny);
  u_ENO1_2d_RK1(:,Ny+ghostcell_width+1:end) = ...
    u_ENO1_2d_RK1(:,ghostcell_width+1:2*ghostcell_width);
  u_ENO1_2d_RK1(1:ghostcell_width,:) = ...
    u_ENO1_2d_RK1(Nx+1:ghostcell_width+Nx,:);
  u_ENO1_2d_RK1(Nx+ghostcell_width+1:end,:) = ...
    u_ENO1_2d_RK1(ghostcell_width+1:2*ghostcell_width,:);
  u_ENO2_2d_RK1(:,1:ghostcell_width) = ...
    u_ENO2_2d_RK1(:,Ny+1:ghostcell_width+Ny);
  u_ENO2_2d_RK1(:,Ny+ghostcell_width+1:end) = ...
    u_ENO2_2d_RK1(:,ghostcell_width+1:2*ghostcell_width);
  u_ENO2_2d_RK1(1:ghostcell_width,:) = ...
    u_ENO2_2d_RK1(Nx+1:ghostcell_width+Nx,:);
  u_ENO2_2d_RK1(Nx+ghostcell_width+1:end,:) = ...
    u_ENO2_2d_RK1(ghostcell_width+1:2*ghostcell_width,:);
  u_ENO3_2d_RK1(:,1:ghostcell_width) = ...
    u_ENO3_2d_RK1(:,Ny+1:ghostcell_width+Ny);
  u_ENO3_2d_RK1(:,Ny+ghostcell_width+1:end) = ...
    u_ENO3_2d_RK1(:,ghostcell_width+1:2*ghostcell_width);
  u_ENO3_2d_RK1(1:ghostcell_width,:) = ...
    u_ENO3_2d_RK1(Nx+1:ghostcell_width+Nx,:);
  u_ENO3_2d_RK1(Nx+ghostcell_width+1:end,:) = ...
    u_ENO3_2d_RK1(ghostcell_width+1:2*ghostcell_width,:);
  u_WENO5_2d_RK1(:,1:ghostcell_width) = ...
    u_WENO5_2d_RK1(:,Ny+1:ghostcell_width+Ny);
  u_WENO5_2d_RK1(:,Ny+ghostcell_width+1:end) = ...
    u_WENO5_2d_RK1(:,ghostcell_width+1:2*ghostcell_width);
  u_WENO5_2d_RK1(1:ghostcell_width,:) = ...
    u_WENO5_2d_RK1(Nx+1:ghostcell_width+Nx,:);
  u_WENO5_2d_RK1(Nx+ghostcell_width+1:end,:) = ...
    u_WENO5_2d_RK1(ghostcell_width+1:2*ghostcell_width,:);

  % compute approximations to u_x
  [u_x_ENO1_2d,u_y_ENO1_2d] = ...
    UPWIND_HJ_ENO1_2D(u_ENO1_2d_RK1,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO2_2d,u_y_ENO2_2d] = ...
    UPWIND_HJ_ENO2_2D(u_ENO2_2d_RK1,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO3_2d,u_y_ENO3_2d] = ...
    UPWIND_HJ_ENO3_2D(u_ENO3_2d_RK1,v_x,v_y,ghostcell_width,dX);
  [u_x_WENO5_2d,u_y_WENO5_2d] = ...
    UPWIND_HJ_WENO5_2D(u_WENO5_2d_RK1,v_x,v_y,ghostcell_width,dX);

  % advance solution
  u_ENO1_2d_RK1 = u_ENO1_2d_RK1 ...
                - dt_RK1*v_x.*u_x_ENO1_2d - dt_RK1*v_y.*u_y_ENO1_2d;
  u_ENO2_2d_RK1 = u_ENO2_2d_RK1 ...
                - dt_RK1*v_x.*u_x_ENO2_2d - dt_RK1*v_y.*u_y_ENO2_2d;
  u_ENO3_2d_RK1 = u_ENO3_2d_RK1 ...
                - dt_RK1*v_x.*u_x_ENO3_2d - dt_RK1*v_y.*u_y_ENO3_2d;
  u_WENO5_2d_RK1 = u_WENO5_2d_RK1 ...
                 - dt_RK1*v_x.*u_x_WENO5_2d - dt_RK1*v_y.*u_y_WENO5_2d;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for TVD RK2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t_RK2)

  % fill boundary cells
  u_ENO1_2d_RK2(:,1:ghostcell_width) = ...
    u_ENO1_2d_RK2(:,Ny+1:ghostcell_width+Ny);
  u_ENO1_2d_RK2(:,Ny+ghostcell_width+1:end) = ...
    u_ENO1_2d_RK2(:,ghostcell_width+1:2*ghostcell_width);
  u_ENO1_2d_RK2(1:ghostcell_width,:) = ...
    u_ENO1_2d_RK2(Nx+1:ghostcell_width+Nx,:);
  u_ENO1_2d_RK2(Nx+ghostcell_width+1:end,:) = ...
    u_ENO1_2d_RK2(ghostcell_width+1:2*ghostcell_width,:);
  u_ENO2_2d_RK2(:,1:ghostcell_width) = ...
    u_ENO2_2d_RK2(:,Ny+1:ghostcell_width+Ny);
  u_ENO2_2d_RK2(:,Ny+ghostcell_width+1:end) = ...
    u_ENO2_2d_RK2(:,ghostcell_width+1:2*ghostcell_width);
  u_ENO2_2d_RK2(1:ghostcell_width,:) = ...
    u_ENO2_2d_RK2(Nx+1:ghostcell_width+Nx,:);
  u_ENO2_2d_RK2(Nx+ghostcell_width+1:end,:) = ...
    u_ENO2_2d_RK2(ghostcell_width+1:2*ghostcell_width,:);
  u_ENO3_2d_RK2(:,1:ghostcell_width) = ...
    u_ENO3_2d_RK2(:,Ny+1:ghostcell_width+Ny);
  u_ENO3_2d_RK2(:,Ny+ghostcell_width+1:end) = ...
    u_ENO3_2d_RK2(:,ghostcell_width+1:2*ghostcell_width);
  u_ENO3_2d_RK2(1:ghostcell_width,:) = ...
    u_ENO3_2d_RK2(Nx+1:ghostcell_width+Nx,:);
  u_ENO3_2d_RK2(Nx+ghostcell_width+1:end,:) = ...
    u_ENO3_2d_RK2(ghostcell_width+1:2*ghostcell_width,:);
  u_WENO5_2d_RK2(:,1:ghostcell_width) = ...
    u_WENO5_2d_RK2(:,Ny+1:ghostcell_width+Ny);
  u_WENO5_2d_RK2(:,Ny+ghostcell_width+1:end) = ...
    u_WENO5_2d_RK2(:,ghostcell_width+1:2*ghostcell_width);
  u_WENO5_2d_RK2(1:ghostcell_width,:) = ...
    u_WENO5_2d_RK2(Nx+1:ghostcell_width+Nx,:);
  u_WENO5_2d_RK2(Nx+ghostcell_width+1:end,:) = ...
    u_WENO5_2d_RK2(ghostcell_width+1:2*ghostcell_width,:);

  % compute approximations to u_x for first stage of TVD RK2
  [u_x_ENO1_2d,u_y_ENO1_2d] = ...
    UPWIND_HJ_ENO1_2D(u_ENO1_2d_RK2,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO2_2d,u_y_ENO2_2d] = ...
    UPWIND_HJ_ENO2_2D(u_ENO2_2d_RK2,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO3_2d,u_y_ENO3_2d] = ...
    UPWIND_HJ_ENO3_2D(u_ENO3_2d_RK2,v_x,v_y,ghostcell_width,dX);
  [u_x_WENO5_2d,u_y_WENO5_2d] = ...
    UPWIND_HJ_WENO5_2D(u_WENO5_2d_RK2,v_x,v_y,ghostcell_width,dX);

  % advance first stage of TVD RK2
  u_ENO1_2d_tmp1 = u_ENO1_2d_RK2 ...
                 - dt_RK2*v_x.*u_x_ENO1_2d - dt_RK2*v_y.*u_y_ENO1_2d;
  u_ENO2_2d_tmp1 = u_ENO2_2d_RK2 ...
                 - dt_RK2*v_x.*u_x_ENO2_2d - dt_RK2*v_y.*u_y_ENO2_2d;
  u_ENO3_2d_tmp1 = u_ENO3_2d_RK2 ...
                 - dt_RK2*v_x.*u_x_ENO3_2d - dt_RK2*v_y.*u_y_ENO3_2d;
  u_WENO5_2d_tmp1 = u_WENO5_2d_RK2 ...
                  - dt_RK2*v_x.*u_x_WENO5_2d - dt_RK2*v_y.*u_y_WENO5_2d;

  % compute approximations to u_x for second stage of TVD RK2
  [u_x_ENO1_2d, u_y_ENO1_2d] = ...
    UPWIND_HJ_ENO1_2D(u_ENO1_2d_tmp1,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO2_2d, u_y_ENO2_2d] = ...
    UPWIND_HJ_ENO2_2D(u_ENO2_2d_tmp1,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO3_2d, u_y_ENO3_2d] = ...
    UPWIND_HJ_ENO3_2D(u_ENO3_2d_tmp1,v_x,v_y,ghostcell_width,dX);
  [u_x_WENO5_2d, u_y_WENO5_2d] = ...
    UPWIND_HJ_WENO5_2D(u_WENO5_2d_tmp1,v_x,v_y,ghostcell_width,dX);

  % advance second stage of TVD RK2
  u_ENO1_2d_RK2 = 0.5*( u_ENO1_2d_RK2 + u_ENO1_2d_tmp1 ...
                      - dt_RK2*v_x.*u_x_ENO1_2d - dt_RK2*v_y.*u_y_ENO1_2d);
  u_ENO2_2d_RK2 = 0.5*( u_ENO2_2d_RK2 + u_ENO2_2d_tmp1 ...
                      - dt_RK2*v_x.*u_x_ENO2_2d - dt_RK2*v_y.*u_y_ENO2_2d);
  u_ENO3_2d_RK2 = 0.5*( u_ENO3_2d_RK2 + u_ENO3_2d_tmp1 ...
                      - dt_RK2*v_x.*u_x_ENO3_2d - dt_RK2*v_y.*u_y_ENO3_2d);
  u_WENO5_2d_RK2 = 0.5*( u_WENO5_2d_RK2 + u_WENO5_2d_tmp1 ...
                       - dt_RK2*v_x.*u_x_WENO5_2d - dt_RK2*v_y.*u_y_WENO5_2d);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for TVD RK3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t_RK3)

  % fill boundary cells
  u_ENO1_2d_RK3(:,1:ghostcell_width) = ...
    u_ENO1_2d_RK3(:,Ny+1:ghostcell_width+Ny);
  u_ENO1_2d_RK3(:,Ny+ghostcell_width+1:end) = ...
    u_ENO1_2d_RK3(:,ghostcell_width+1:2*ghostcell_width);
  u_ENO1_2d_RK3(1:ghostcell_width,:) = ...
    u_ENO1_2d_RK3(Nx+1:ghostcell_width+Nx,:);
  u_ENO1_2d_RK3(Nx+ghostcell_width+1:end,:) = ...
    u_ENO1_2d_RK3(ghostcell_width+1:2*ghostcell_width,:);
  u_ENO2_2d_RK3(:,1:ghostcell_width) = ...
    u_ENO2_2d_RK3(:,Ny+1:ghostcell_width+Ny);
  u_ENO2_2d_RK3(:,Ny+ghostcell_width+1:end) = ...
    u_ENO2_2d_RK3(:,ghostcell_width+1:2*ghostcell_width);
  u_ENO2_2d_RK3(1:ghostcell_width,:) = ...
    u_ENO2_2d_RK3(Nx+1:ghostcell_width+Nx,:);
  u_ENO2_2d_RK3(Nx+ghostcell_width+1:end,:) = ...
    u_ENO2_2d_RK3(ghostcell_width+1:2*ghostcell_width,:);
  u_ENO3_2d_RK3(:,1:ghostcell_width) = ...
    u_ENO3_2d_RK3(:,Ny+1:ghostcell_width+Ny);
  u_ENO3_2d_RK3(:,Ny+ghostcell_width+1:end) = ...
    u_ENO3_2d_RK3(:,ghostcell_width+1:2*ghostcell_width);
  u_ENO3_2d_RK3(1:ghostcell_width,:) = ...
    u_ENO3_2d_RK3(Nx+1:ghostcell_width+Nx,:);
  u_ENO3_2d_RK3(Nx+ghostcell_width+1:end,:) = ...
    u_ENO3_2d_RK3(ghostcell_width+1:2*ghostcell_width,:);
  u_WENO5_2d_RK3(:,1:ghostcell_width) = ...
    u_WENO5_2d_RK3(:,Ny+1:ghostcell_width+Ny);
  u_WENO5_2d_RK3(:,Ny+ghostcell_width+1:end) = ...
    u_WENO5_2d_RK3(:,ghostcell_width+1:2*ghostcell_width);
  u_WENO5_2d_RK3(1:ghostcell_width,:) = ...
    u_WENO5_2d_RK3(Nx+1:ghostcell_width+Nx,:);
  u_WENO5_2d_RK3(Nx+ghostcell_width+1:end,:) = ...
    u_WENO5_2d_RK3(ghostcell_width+1:2*ghostcell_width,:);

  % compute approximations to u_x
  [u_x_ENO1_2d, u_y_ENO1_2d] = ...
    UPWIND_HJ_ENO1_2D(u_ENO1_2d_RK3,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO2_2d, u_y_ENO2_2d] = ...
    UPWIND_HJ_ENO2_2D(u_ENO2_2d_RK3,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO3_2d, u_y_ENO3_2d] = ...
    UPWIND_HJ_ENO3_2D(u_ENO3_2d_RK3,v_x,v_y,ghostcell_width,dX);
  [u_x_WENO5_2d, u_y_WENO5_2d] = ...
    UPWIND_HJ_WENO5_2D(u_WENO5_2d_RK3,v_x,v_y,ghostcell_width,dX);

  % advance first stage of TVD RK3
  u_ENO1_2d_tmp1 = u_ENO1_2d_RK3 ...
                 - dt_RK3*v_x.*u_x_ENO1_2d - dt_RK3*v_y.*u_y_ENO1_2d;
  u_ENO2_2d_tmp1 = u_ENO2_2d_RK3 ...
                 - dt_RK3*v_x.*u_x_ENO2_2d - dt_RK3*v_y.*u_y_ENO2_2d;
  u_ENO3_2d_tmp1 = u_ENO3_2d_RK3 ...
                 - dt_RK3*v_x.*u_x_ENO3_2d - dt_RK3*v_y.*u_y_ENO3_2d;
  u_WENO5_2d_tmp1 = u_WENO5_2d_RK3 ...
                  - dt_RK3*v_x.*u_x_WENO5_2d - dt_RK3*v_y.*u_y_WENO5_2d;

  % compute approximations to u_x for second stage of TVD RK3
  [u_x_ENO1_2d, u_y_ENO1_2d] = ...
    UPWIND_HJ_ENO1_2D(u_ENO1_2d_tmp1,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO2_2d, u_y_ENO2_2d] = ...
    UPWIND_HJ_ENO2_2D(u_ENO2_2d_tmp1,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO3_2d, u_y_ENO3_2d] = ...
    UPWIND_HJ_ENO3_2D(u_ENO3_2d_tmp1,v_x,v_y,ghostcell_width,dX);
  [u_x_WENO5_2d, u_y_WENO5_2d] = ...
    UPWIND_HJ_WENO5_2D(u_WENO5_2d_tmp1,v_x,v_y,ghostcell_width,dX);

  % advance second stage of TVD RK3
  u_ENO1_2d_tmp2 = 0.75*u_ENO1_2d_RK3 ...
                 + 0.25*( u_ENO1_2d_tmp1 ...
                        - dt_RK3*v_x.*u_x_ENO1_2d ...
                        - dt_RK3*v_y.*u_y_ENO1_2d );
  u_ENO2_2d_tmp2 = 0.75*u_ENO2_2d_RK3 ...
                 + 0.25*( u_ENO2_2d_tmp1 ...
                        - dt_RK3*v_x.*u_x_ENO2_2d ...
                        - dt_RK3*v_y.*u_y_ENO2_2d );
  u_ENO3_2d_tmp2 = 0.75*u_ENO3_2d_RK3 ...
                 + 0.25*( u_ENO3_2d_tmp1 ...
                        - dt_RK3*v_x.*u_x_ENO3_2d ...
                        - dt_RK3*v_y.*u_y_ENO3_2d );
  u_WENO5_2d_tmp2 = 0.75*u_WENO5_2d_RK3 ...
                  + 0.25*( u_WENO5_2d_tmp1 ...
                         - dt_RK3*v_x.*u_x_WENO5_2d ...
                         - dt_RK3*v_y.*u_y_WENO5_2d );

  % compute approximations to u_x for third stage of TVD RK3
  [u_x_ENO1_2d, u_y_ENO1_2d] = ...
    UPWIND_HJ_ENO1_2D(u_ENO1_2d_tmp2,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO2_2d, u_y_ENO2_2d] = ...
    UPWIND_HJ_ENO2_2D(u_ENO2_2d_tmp2,v_x,v_y,ghostcell_width,dX);
  [u_x_ENO3_2d, u_y_ENO3_2d] = ...
    UPWIND_HJ_ENO3_2D(u_ENO3_2d_tmp2,v_x,v_y,ghostcell_width,dX);
  [u_x_WENO5_2d, u_y_WENO5_2d] = ...
    UPWIND_HJ_WENO5_2D(u_WENO5_2d_tmp2,v_x,v_y,ghostcell_width,dX);

  % advance third stage of TVD RK3
  u_ENO1_2d_RK3 = 1/3*u_ENO1_2d_RK3 ...
                + 2/3*( u_ENO1_2d_tmp2 ...
                      - dt_RK3*v_x.*u_x_ENO1_2d ...
                      - dt_RK3*v_y.*u_y_ENO1_2d);
  u_ENO2_2d_RK3 = 1/3*u_ENO2_2d_RK3 ...
                + 2/3*( u_ENO2_2d_tmp2 ...
                      - dt_RK3*v_x.*u_x_ENO2_2d ...
                      - dt_RK3*v_y.*u_y_ENO2_2d);
  u_ENO3_2d_RK3 = 1/3*u_ENO3_2d_RK3 ...
                + 2/3*( u_ENO3_2d_tmp2 ...
                      - dt_RK3*v_x.*u_x_ENO3_2d ...
                      - dt_RK3*v_y.*u_y_ENO3_2d);
  u_WENO5_2d_RK3 = 1/3*u_WENO5_2d_RK3 ...
                 + 2/3*( u_WENO5_2d_tmp2 ...
                       - dt_RK3*v_x.*u_x_WENO5_2d ...
                       - dt_RK3*v_y.*u_y_WENO5_2d);

end

% plot results
figure(1); clf;
surf(x,y,u_ENO3_2d_RK2);
xlabel('x');
ylabel('y');
axis([-1 1 -1 1]);

return
plot(x,u_ENO2_2d_RK1,'g-.');
plot(x,u_ENO3_2d_RK1,'r-.');
plot(x,u_WENO5_2d_RK1,'m-.');

plot(x,u_ENO1_2d_RK2,'b--');
plot(x,u_ENO2_2d_RK2,'g--');
plot(x,u_ENO3_2d_RK2,'r--');
plot(x,u_WENO5_2d_RK2,'m--');

plot(x,u_ENO1_2d_RK3,'b');
plot(x,u_ENO2_2d_RK3,'g');
plot(x,u_ENO3_2d_RK3,'r');
plot(x,u_WENO5_2d_RK3,'m');

