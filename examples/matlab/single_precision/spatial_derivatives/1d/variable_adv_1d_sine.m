%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File:        variable_adv_1d_sine.m
% Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:    $Revision: 1.3 $
% Modified:    $Date: 2006/01/24 21:45:49 $
% Description: MATLAB test program for 1D ENO/WENO spatial derivatives 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script tests the 1D ENO/WENO derivative calculations using 
% the advection equation:
%
%    u_t + v(x) u_x = 0
%
% where v(x) is a position dependent advection velocity.  The initial
% condition is a sine wave, u(x,t=0) = sin(pi*x) and the boundary
% conditions are periodic.
%
% In this code, time advection is done using forward euler (TVD RK1), 
% TVD RK2, and TVD RK3.
%
% Kevin Chu
% Dept of Mathematics, MIT
% March 2005
%

% setup environment
clear
format long

% set up spatial grid parameters
N = 500;
max_ENO_order = 3;
ghostcell_width = max_ENO_order;
N_with_ghostcells = N+2*ghostcell_width;
x_lo = -1;
x_hi = 1;
dx = (x_hi-x_lo)/N;
x = (x_lo-(ghostcell_width-0.5)*dx:dx:x_hi+ghostcell_width*dx)';

% set advection velocity function
vel = single( (sin(pi*x)+1.5) );

% set up time integration parameters
cfl_RK1 = 0.3;
cfl_RK2 = 0.7;
cfl_RK3 = 1;
t_i = 0;
t_f = 5;
dt_RK1 = cfl_RK1*dx/norm(vel,inf);
dt_RK2 = cfl_RK2*dx/norm(vel,inf);
	dt_RK3 = cfl_RK3*dx/norm(vel,inf);
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
	u_ENO1_1d_RK1 = single( sin(2*pi*x) );
u_ENO2_1d_RK1 = u_ENO1_1d_RK1;
u_ENO3_1d_RK1 = u_ENO1_1d_RK1;
u_WENO5_1d_RK1 = u_ENO1_1d_RK1;
u_ENO1_1d_RK2 = u_ENO1_1d_RK1;
u_ENO2_1d_RK2 = u_ENO1_1d_RK1;
u_ENO3_1d_RK2 = u_ENO1_1d_RK1;
u_WENO5_1d_RK2 = u_ENO1_1d_RK1;
u_ENO1_1d_RK3 = u_ENO1_1d_RK1;
u_ENO2_1d_RK3 = u_ENO1_1d_RK1;
u_ENO3_1d_RK3 = u_ENO1_1d_RK1;
u_WENO5_1d_RK3 = u_ENO1_1d_RK1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for forward euler (TVD RK1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t_RK1)

  % fill boundary cells
  u_ENO1_1d_RK1(1:ghostcell_width) = u_ENO1_1d_RK1(N+1:ghostcell_width+N);
  u_ENO1_1d_RK1(N+ghostcell_width+1:end) = ...
    u_ENO1_1d_RK1(ghostcell_width+1:2*ghostcell_width);
  u_ENO2_1d_RK1(1:ghostcell_width) = u_ENO2_1d_RK1(N+1:ghostcell_width+N);
  u_ENO2_1d_RK1(N+ghostcell_width+1:end) = ...
    u_ENO2_1d_RK1(ghostcell_width+1:2*ghostcell_width);
  u_ENO3_1d_RK1(1:ghostcell_width) = u_ENO3_1d_RK1(N+1:ghostcell_width+N);
  u_ENO3_1d_RK1(N+ghostcell_width+1:end) = ...
    u_ENO3_1d_RK1(ghostcell_width+1:2*ghostcell_width);
  u_WENO5_1d_RK1(1:ghostcell_width) = u_WENO5_1d_RK1(N+1:ghostcell_width+N);
  u_WENO5_1d_RK1(N+ghostcell_width+1:end) = ...
    u_WENO5_1d_RK1(ghostcell_width+1:2*ghostcell_width);

  % compute approximations to u_x
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_RK1,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_RK1,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_RK1,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_RK1,vel,ghostcell_width,dx);

  % advance solution
  u_ENO1_1d_RK1 = u_ENO1_1d_RK1 - dt_RK1*vel.*u_x_ENO1_1d;
  u_ENO2_1d_RK1 = u_ENO2_1d_RK1 - dt_RK1*vel.*u_x_ENO2_1d;
  u_ENO3_1d_RK1 = u_ENO3_1d_RK1 - dt_RK1*vel.*u_x_ENO3_1d;
  u_WENO5_1d_RK1 = u_WENO5_1d_RK1 - dt_RK1*vel.*u_x_WENO5_1d;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for TVD RK2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t_RK2)

  % fill boundary cells
  u_ENO1_1d_RK2(1:ghostcell_width) = u_ENO1_1d_RK2(N+1:ghostcell_width+N);
  u_ENO1_1d_RK2(N+ghostcell_width+1:end) = ...
    u_ENO1_1d_RK2(ghostcell_width+1:2*ghostcell_width);
  u_ENO2_1d_RK2(1:ghostcell_width) = u_ENO2_1d_RK2(N+1:ghostcell_width+N);
  u_ENO2_1d_RK2(N+ghostcell_width+1:end) = ...
    u_ENO2_1d_RK2(ghostcell_width+1:2*ghostcell_width);
  u_ENO3_1d_RK2(1:ghostcell_width) = u_ENO3_1d_RK2(N+1:ghostcell_width+N);
  u_ENO3_1d_RK2(N+ghostcell_width+1:end) = ...
    u_ENO3_1d_RK2(ghostcell_width+1:2*ghostcell_width);
  u_WENO5_1d_RK2(1:ghostcell_width) = u_WENO5_1d_RK2(N+1:ghostcell_width+N);
  u_WENO5_1d_RK2(N+ghostcell_width+1:end) = ...
    u_WENO5_1d_RK2(ghostcell_width+1:2*ghostcell_width);

  % compute approximations to u_x for first stage of TVD RK2
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_RK2,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_RK2,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_RK2,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_RK2,vel,ghostcell_width,dx);

  % advance first stage of TVD RK2
  u_ENO1_1d_tmp1 = u_ENO1_1d_RK2 - dt_RK2*vel.*u_x_ENO1_1d;
  u_ENO2_1d_tmp1 = u_ENO2_1d_RK2 - dt_RK2*vel.*u_x_ENO2_1d;
  u_ENO3_1d_tmp1 = u_ENO3_1d_RK2 - dt_RK2*vel.*u_x_ENO3_1d;
  u_WENO5_1d_tmp1 = u_WENO5_1d_RK2 - dt_RK2*vel.*u_x_WENO5_1d;

  % compute approximations to u_x for second stage of TVD RK2
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_tmp1,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_tmp1,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_tmp1,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_tmp1,vel,ghostcell_width,dx);

  % advance second stage of TVD RK2
  u_ENO1_1d_RK2 = 0.5*(u_ENO1_1d_RK2 + u_ENO1_1d_tmp1 - dt_RK2*vel.*u_x_ENO1_1d);
  u_ENO2_1d_RK2 = 0.5*(u_ENO2_1d_RK2 + u_ENO2_1d_tmp1 - dt_RK2*vel.*u_x_ENO2_1d);
  u_ENO3_1d_RK2 = 0.5*(u_ENO3_1d_RK2 + u_ENO3_1d_tmp1 - dt_RK2*vel.*u_x_ENO3_1d);
  u_WENO5_1d_RK2 = 0.5*(u_WENO5_1d_RK2 + u_WENO5_1d_tmp1 - dt_RK2*vel.*u_x_WENO5_1d);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for TVD RK3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t_RK3)

  % fill boundary cells
  u_ENO1_1d_RK3(1:ghostcell_width) = u_ENO1_1d_RK3(N+1:ghostcell_width+N);
  u_ENO1_1d_RK3(N+ghostcell_width+1:end) = ...
    u_ENO1_1d_RK3(ghostcell_width+1:2*ghostcell_width);
  u_ENO2_1d_RK3(1:ghostcell_width) = u_ENO2_1d_RK3(N+1:ghostcell_width+N);
  u_ENO2_1d_RK3(N+ghostcell_width+1:end) = ...
    u_ENO2_1d_RK3(ghostcell_width+1:2*ghostcell_width);
  u_ENO3_1d_RK3(1:ghostcell_width) = u_ENO3_1d_RK3(N+1:ghostcell_width+N);
  u_ENO3_1d_RK3(N+ghostcell_width+1:end) = ...
    u_ENO3_1d_RK3(ghostcell_width+1:2*ghostcell_width);
  u_WENO5_1d_RK3(1:ghostcell_width) = u_WENO5_1d_RK3(N+1:ghostcell_width+N);
  u_WENO5_1d_RK3(N+ghostcell_width+1:end) = ...
    u_WENO5_1d_RK3(ghostcell_width+1:2*ghostcell_width);

  % compute approximations to u_x
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_RK3,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_RK3,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_RK3,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_RK3,vel,ghostcell_width,dx);

  % advance first stage of TVD RK3
  u_ENO1_1d_tmp1 = u_ENO1_1d_RK3 - dt_RK3*vel.*u_x_ENO1_1d;
  u_ENO2_1d_tmp1 = u_ENO2_1d_RK3 - dt_RK3*vel.*u_x_ENO2_1d;
  u_ENO3_1d_tmp1 = u_ENO3_1d_RK3 - dt_RK3*vel.*u_x_ENO3_1d;
  u_WENO5_1d_tmp1 = u_WENO5_1d_RK3 - dt_RK3*vel.*u_x_WENO5_1d;

  % compute approximations to u_x for second stage of TVD RK3
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_tmp1,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_tmp1,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_tmp1,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_tmp1,vel,ghostcell_width,dx);

  % advance second stage of TVD RK3
  u_ENO1_1d_tmp2 = 0.75*u_ENO1_1d_RK3 + 0.25*(u_ENO1_1d_tmp1 - dt_RK3*vel.*u_x_ENO1_1d);
  u_ENO2_1d_tmp2 = 0.75*u_ENO2_1d_RK3 + 0.25*(u_ENO2_1d_tmp1 - dt_RK3*vel.*u_x_ENO2_1d);
  u_ENO3_1d_tmp2 = 0.75*u_ENO3_1d_RK3 + 0.25*(u_ENO3_1d_tmp1 - dt_RK3*vel.*u_x_ENO3_1d);
  u_WENO5_1d_tmp2 = 0.75*u_WENO5_1d_RK3 + 0.25*(u_WENO5_1d_tmp1 - dt_RK3*vel.*u_x_WENO5_1d);

  % compute approximations to u_x for third stage of TVD RK3
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_tmp2,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_tmp2,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_tmp2,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_tmp2,vel,ghostcell_width,dx);

  % advance third stage of TVD RK3
  u_ENO1_1d_RK3 = 1/3*u_ENO1_1d_RK3 + 2/3*(u_ENO1_1d_tmp2 - dt_RK3*vel.*u_x_ENO1_1d);
  u_ENO2_1d_RK3 = 1/3*u_ENO2_1d_RK3 + 2/3*(u_ENO2_1d_tmp2 - dt_RK3*vel.*u_x_ENO2_1d);
  u_ENO3_1d_RK3 = 1/3*u_ENO3_1d_RK3 + 2/3*(u_ENO3_1d_tmp2 - dt_RK3*vel.*u_x_ENO3_1d);
  u_WENO5_1d_RK3 = 1/3*u_WENO5_1d_RK3 + 2/3*(u_WENO5_1d_tmp2 - dt_RK3*vel.*u_x_WENO5_1d);

end

% plot results
figure(1); clf;
plot(x,u_ENO1_1d_RK1,'b-.');
hold on;
plot(x,u_ENO2_1d_RK1,'g-.');
plot(x,u_ENO3_1d_RK1,'r-.');
plot(x,u_WENO5_1d_RK1,'m-.');

plot(x,u_ENO1_1d_RK2,'b--');
plot(x,u_ENO2_1d_RK2,'g--');
plot(x,u_ENO3_1d_RK2,'r--');
plot(x,u_WENO5_1d_RK2,'m--');

plot(x,u_ENO1_1d_RK3,'b');
plot(x,u_ENO2_1d_RK3,'g');
plot(x,u_ENO3_1d_RK3,'r');
plot(x,u_WENO5_1d_RK3,'m');
axis([-1 1 -1.1 1.1]);

