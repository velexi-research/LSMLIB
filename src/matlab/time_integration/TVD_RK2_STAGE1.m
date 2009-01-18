%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TVD_RK2_STAGE1 advances the solution through the first stage of a 
%                second-order TVD Runge-Kutta step.
%
% Usage: u_stage1 = TVD_RK2_STAGE1(u_cur, rhs, dt)
%
% Arguments:
% - u_cur:       u(t_cur)
% - rhs:         right-hand side of time evolution equation at t_cur
% - dt:          step size
%
% Return value:
% - u_stage1:    u_approx(t_cur + dt)
%
% NOTES:
% - u_cur and rhs are assumed to have the same dimensions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright:  (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:   $Revision$
% Modified:   $Date$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_stage1 = TVD_RK2_STAGE1(u_cur, rhs, dt)

u_stage1 = u_cur + dt*rhs;
