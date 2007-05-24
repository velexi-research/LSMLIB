%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TVD_RK1_STEP takes a single first-order Runge-Kutta (i.e.
%              Forward Euler) step.
%
% Usage: u_next = TVD_RK1_STEP(u_cur, rhs, dt)
%
% Arguments:
% - u_cur:       u(t_cur)
% - rhs:         right-hand side of time evolution equation at t_cur
% - dt:          step size
%
% Return value:
% - u_next:      u(t_cur + dt)
%
% NOTES:
% - u_cur and rhs are assumed to have the same dimensions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author:     Kevin T. Chu 
% Copyright:  (c) 2005-2006, MAE Princeton University 
% Revision:   $Revision: 1.3 $
% Modified:   $Date: 2006/04/22 12:39:34 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_next = TVD_RK1_STEP(u_cur, rhs, dt)

u_next = u_cur + dt*rhs;
