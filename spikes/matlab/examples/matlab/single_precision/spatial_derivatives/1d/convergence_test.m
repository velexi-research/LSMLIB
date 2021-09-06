%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File:        convergence_test.m
% Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
%                  Regents of the University of Texas.  All rights reserved.
%              (c) 2009 Kevin T. Chu.  All rights reserved.
% Revision:    $Revision$
% Modified:    $Date$
% Description: MATLAB test program for 1D ENO/WENO spatial derivatives 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script tests the 1D ENO derivative calculations using the
% advection equation:
%
%    u_t + v(x) u_x = 0
%
% where v(x) is a position dependent advection velocity.  The initial
% condition is a sine wave, u(x,t=0) = sin(pi*x) and the boundary
% conditions are periodic.
%
% In this code, time advection is done using forward euler (TVD RK1), 
% TVD RK2, and TVD RK3_low_res.
%
% Kevin Chu
% Dept of Mathematics, MIT
% March 2005
%

% setup environment
clear
format long

% set up spatial grid parameters
N = 200;
max_ENO_order = 3;
ghostcell_width = max_ENO_order;
N_with_ghostcells = N+2*ghostcell_width;
x_lo = -1;
x_hi = 1;
dx = (x_hi-x_lo)/N;
x_low_res = [x_lo-(ghostcell_width-0.5)*dx:dx:x_hi+ghostcell_width*dx]';

% set advection velocity function
vel = single( (sin(pi*(x_low_res-0.25))) );
vel = single( -(sin(pi*x_low_res)+2) );

% set up time integration parameters
cfl_RK3_low_res = 0.5/norm(vel,inf);
t_i = 0;
t_f = 1;
dt_RK3_low_res = cfl_RK3_low_res*dx;
t_RK3_low_res = t_i:dt_RK3_low_res:t_f;
if (t_RK3_low_res(end) ~= t_f)
  t_RK3_low_res = [t_RK3_low_res t_f];
end

% initialize u
u_ENO1_RK3_low_res = single( sin(2*pi*x_low_res) );
u_ENO2_RK3_low_res = u_ENO1_RK3_low_res;
u_ENO3_RK3_low_res = u_ENO1_RK3_low_res;
u_WENO5_RK3_low_res = u_ENO1_RK3_low_res;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for TVD RK3_low_res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t_RK3_low_res)

  % fill boundary cells
  u_ENO1_RK3_low_res(1:ghostcell_width) = u_ENO1_RK3_low_res(N+1:ghostcell_width+N);
  u_ENO1_RK3_low_res(N+ghostcell_width+1:end) = ...
    u_ENO1_RK3_low_res(ghostcell_width+1:2*ghostcell_width);
  u_ENO2_RK3_low_res(1:ghostcell_width) = u_ENO2_RK3_low_res(N+1:ghostcell_width+N);
  u_ENO2_RK3_low_res(N+ghostcell_width+1:end) = ...
    u_ENO2_RK3_low_res(ghostcell_width+1:2*ghostcell_width);
  u_ENO3_RK3_low_res(1:ghostcell_width) = u_ENO3_RK3_low_res(N+1:ghostcell_width+N);
  u_ENO3_RK3_low_res(N+ghostcell_width+1:end) = ...
    u_ENO3_RK3_low_res(ghostcell_width+1:2*ghostcell_width);
  u_WENO5_RK3_low_res(1:ghostcell_width) = u_WENO5_RK3_low_res(N+1:ghostcell_width+N);
  u_WENO5_RK3_low_res(N+ghostcell_width+1:end) = ...
    u_WENO5_RK3_low_res(ghostcell_width+1:2*ghostcell_width);

  % compute approximations to u_x
  u_x_ENO1 = UPWIND_HJ_ENO1_1D(u_ENO1_RK3_low_res,vel,ghostcell_width,dx);
  u_x_ENO2 = UPWIND_HJ_ENO2_1D(u_ENO2_RK3_low_res,vel,ghostcell_width,dx);
  u_x_ENO3 = UPWIND_HJ_ENO3_1D(u_ENO3_RK3_low_res,vel,ghostcell_width,dx);
  u_x_WENO5 = UPWIND_HJ_WENO5_1D(u_WENO5_RK3_low_res,vel,ghostcell_width,dx);

  % advance first stage of TVD RK3_low_res
  u_ENO1_tmp1 = u_ENO1_RK3_low_res - dt_RK3_low_res*vel.*u_x_ENO1;
  u_ENO2_tmp1 = u_ENO2_RK3_low_res - dt_RK3_low_res*vel.*u_x_ENO2;
  u_ENO3_tmp1 = u_ENO3_RK3_low_res - dt_RK3_low_res*vel.*u_x_ENO3;
  u_WENO5_tmp1 = u_WENO5_RK3_low_res - dt_RK3_low_res*vel.*u_x_WENO5;

  % compute approximations to u_x for second stage of TVD RK3_low_res
  u_x_ENO1 = UPWIND_HJ_ENO1_1D(u_ENO1_tmp1,vel,ghostcell_width,dx);
  u_x_ENO2 = UPWIND_HJ_ENO2_1D(u_ENO2_tmp1,vel,ghostcell_width,dx);
  u_x_ENO3 = UPWIND_HJ_ENO3_1D(u_ENO3_tmp1,vel,ghostcell_width,dx);
  u_x_WENO5 = UPWIND_HJ_WENO5_1D(u_WENO5_tmp1,vel,ghostcell_width,dx);

  % advance second stage of TVD RK3_low_res
  u_ENO1_tmp2 = 0.75*u_ENO1_RK3_low_res + 0.25*(u_ENO1_tmp1 - dt_RK3_low_res*vel.*u_x_ENO1);
  u_ENO2_tmp2 = 0.75*u_ENO2_RK3_low_res + 0.25*(u_ENO2_tmp1 - dt_RK3_low_res*vel.*u_x_ENO2);
  u_ENO3_tmp2 = 0.75*u_ENO3_RK3_low_res + 0.25*(u_ENO3_tmp1 - dt_RK3_low_res*vel.*u_x_ENO3);
  u_WENO5_tmp2 = 0.75*u_WENO5_RK3_low_res + 0.25*(u_WENO5_tmp1 - dt_RK3_low_res*vel.*u_x_WENO5);

  % compute approximations to u_x for third stage of TVD RK3_low_res
  u_x_ENO1 = UPWIND_HJ_ENO1_1D(u_ENO1_tmp2,vel,ghostcell_width,dx);
  u_x_ENO2 = UPWIND_HJ_ENO2_1D(u_ENO2_tmp2,vel,ghostcell_width,dx);
  u_x_ENO3 = UPWIND_HJ_ENO3_1D(u_ENO3_tmp2,vel,ghostcell_width,dx);
  u_x_WENO5 = UPWIND_HJ_WENO5_1D(u_WENO5_tmp2,vel,ghostcell_width,dx);

  % advance third stage of TVD RK3_low_res
  u_ENO1_RK3_low_res = 1/3*u_ENO1_RK3_low_res + 2/3*(u_ENO1_tmp2 - dt_RK3_low_res*vel.*u_x_ENO1);
  u_ENO2_RK3_low_res = 1/3*u_ENO2_RK3_low_res + 2/3*(u_ENO2_tmp2 - dt_RK3_low_res*vel.*u_x_ENO2);
  u_ENO3_RK3_low_res = 1/3*u_ENO3_RK3_low_res + 2/3*(u_ENO3_tmp2 - dt_RK3_low_res*vel.*u_x_ENO3);
  u_WENO5_RK3_low_res = 1/3*u_WENO5_RK3_low_res + 2/3*(u_WENO5_tmp2 - dt_RK3_low_res*vel.*u_x_WENO5);

end


% set up spatial grid parameters
N = 500;
max_ENO_order = 3;
ghostcell_width = max_ENO_order;
N_with_ghostcells = N+2*ghostcell_width;
x_lo = -1;
x_hi = 1;
dx = (x_hi-x_lo)/N;
x_hi_res = [x_lo-(ghostcell_width-0.5)*dx:dx:x_hi+ghostcell_width*dx]';

% set advection velocity function
vel = single( (sin(pi*(x_hi_res-0.25))) );
vel = single( -(sin(pi*x_hi_res)+2) );

% set up time integration parameters
cfl_RK3_hi_res = 0.5/norm(vel,inf);
t_i = 0;
t_f = 1;
dt_RK3_hi_res = cfl_RK3_hi_res*dx;
t_RK3_hi_res = t_i:dt_RK3_hi_res:t_f;
if (t_RK3_hi_res(end) ~= t_f)
  t_RK3_hi_res = [t_RK3_hi_res t_f];
end

% initialize u
u_ENO1_RK3_hi_res = single( sin(2*pi*x_hi_res) );
u_ENO2_RK3_hi_res = u_ENO1_RK3_hi_res;
u_ENO3_RK3_hi_res = u_ENO1_RK3_hi_res;
u_WENO5_RK3_hi_res = u_ENO1_RK3_hi_res;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for TVD RK3_hi_res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t_RK3_hi_res)

  % fill boundary cells
  u_ENO1_RK3_hi_res(1:ghostcell_width) = u_ENO1_RK3_hi_res(N+1:ghostcell_width+N);
  u_ENO1_RK3_hi_res(N+ghostcell_width+1:end) = ...
    u_ENO1_RK3_hi_res(ghostcell_width+1:2*ghostcell_width);
  u_ENO2_RK3_hi_res(1:ghostcell_width) = u_ENO2_RK3_hi_res(N+1:ghostcell_width+N);
  u_ENO2_RK3_hi_res(N+ghostcell_width+1:end) = ...
    u_ENO2_RK3_hi_res(ghostcell_width+1:2*ghostcell_width);
  u_ENO3_RK3_hi_res(1:ghostcell_width) = u_ENO3_RK3_hi_res(N+1:ghostcell_width+N);
  u_ENO3_RK3_hi_res(N+ghostcell_width+1:end) = ...
    u_ENO3_RK3_hi_res(ghostcell_width+1:2*ghostcell_width);
  u_WENO5_RK3_hi_res(1:ghostcell_width) = u_WENO5_RK3_hi_res(N+1:ghostcell_width+N);
  u_WENO5_RK3_hi_res(N+ghostcell_width+1:end) = ...
    u_WENO5_RK3_hi_res(ghostcell_width+1:2*ghostcell_width);

  % compute approximations to u_x
  u_x_ENO1 = UPWIND_HJ_ENO1_1D(u_ENO1_RK3_hi_res,vel,ghostcell_width,dx);
  u_x_ENO2 = UPWIND_HJ_ENO2_1D(u_ENO2_RK3_hi_res,vel,ghostcell_width,dx);
  u_x_ENO3 = UPWIND_HJ_ENO3_1D(u_ENO3_RK3_hi_res,vel,ghostcell_width,dx);
  u_x_WENO5 = UPWIND_HJ_WENO5_1D(u_WENO5_RK3_hi_res,vel,ghostcell_width,dx);

  % advance first stage of TVD RK3_hi_res
  u_ENO1_tmp1 = u_ENO1_RK3_hi_res - dt_RK3_hi_res*vel.*u_x_ENO1;
  u_ENO2_tmp1 = u_ENO2_RK3_hi_res - dt_RK3_hi_res*vel.*u_x_ENO2;
  u_ENO3_tmp1 = u_ENO3_RK3_hi_res - dt_RK3_hi_res*vel.*u_x_ENO3;
  u_WENO5_tmp1 = u_WENO5_RK3_hi_res - dt_RK3_hi_res*vel.*u_x_WENO5;

  % compute approximations to u_x for second stage of TVD RK3_hi_res
  u_x_ENO1 = UPWIND_HJ_ENO1_1D(u_ENO1_tmp1,vel,ghostcell_width,dx);
  u_x_ENO2 = UPWIND_HJ_ENO2_1D(u_ENO2_tmp1,vel,ghostcell_width,dx);
  u_x_ENO3 = UPWIND_HJ_ENO3_1D(u_ENO3_tmp1,vel,ghostcell_width,dx);
  u_x_WENO5 = UPWIND_HJ_WENO5_1D(u_WENO5_tmp1,vel,ghostcell_width,dx);

  % advance second stage of TVD RK3_hi_res
  u_ENO1_tmp2 = 0.75*u_ENO1_RK3_hi_res + 0.25*(u_ENO1_tmp1 - dt_RK3_hi_res*vel.*u_x_ENO1);
  u_ENO2_tmp2 = 0.75*u_ENO2_RK3_hi_res + 0.25*(u_ENO2_tmp1 - dt_RK3_hi_res*vel.*u_x_ENO2);
  u_ENO3_tmp2 = 0.75*u_ENO3_RK3_hi_res + 0.25*(u_ENO3_tmp1 - dt_RK3_hi_res*vel.*u_x_ENO3);
  u_WENO5_tmp2 = 0.75*u_WENO5_RK3_hi_res + 0.25*(u_WENO5_tmp1 - dt_RK3_hi_res*vel.*u_x_WENO5);

  % compute approximations to u_x for third stage of TVD RK3_low_res
  u_x_ENO1 = UPWIND_HJ_ENO1_1D(u_ENO1_tmp2,vel,ghostcell_width,dx);
  u_x_ENO2 = UPWIND_HJ_ENO2_1D(u_ENO2_tmp2,vel,ghostcell_width,dx);
  u_x_ENO3 = UPWIND_HJ_ENO3_1D(u_ENO3_tmp2,vel,ghostcell_width,dx);
  u_x_WENO5 = UPWIND_HJ_WENO5_1D(u_WENO5_tmp2,vel,ghostcell_width,dx);

  % advance third stage of TVD RK3_hi_res
  u_ENO1_RK3_hi_res = 1/3*u_ENO1_RK3_hi_res + 2/3*(u_ENO1_tmp2 - dt_RK3_hi_res*vel.*u_x_ENO1);
  u_ENO2_RK3_hi_res = 1/3*u_ENO2_RK3_hi_res + 2/3*(u_ENO2_tmp2 - dt_RK3_hi_res*vel.*u_x_ENO2);
  u_ENO3_RK3_hi_res = 1/3*u_ENO3_RK3_hi_res + 2/3*(u_ENO3_tmp2 - dt_RK3_hi_res*vel.*u_x_ENO3);
  u_WENO5_RK3_hi_res = 1/3*u_WENO5_RK3_hi_res + 2/3*(u_WENO5_tmp2 - dt_RK3_hi_res*vel.*u_x_WENO5);

end

% plot results
figure(1); clf;
plot(x_low_res,u_ENO1_RK3_low_res,'b--');
hold on;
plot(x_low_res,u_ENO2_RK3_low_res,'g--');
plot(x_low_res,u_ENO3_RK3_low_res,'r--');
plot(x_low_res,u_WENO5_RK3_low_res,'m--');
plot(x_hi_res,u_ENO1_RK3_hi_res,'b-');
plot(x_hi_res,u_ENO2_RK3_hi_res,'g-');
plot(x_hi_res,u_ENO3_RK3_hi_res,'r-');
plot(x_hi_res,u_WENO5_RK3_hi_res,'m-');
axis([-1 1 -1.1 1.1]);
