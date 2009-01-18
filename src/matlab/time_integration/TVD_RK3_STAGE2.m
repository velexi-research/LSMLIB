%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TVD_RK3_STAGE2 advances the solution through the second stage of a 
%                third-order TVD Runge-Kutta step.
%
% Usage: u_stage2 = TVD_RK3_STAGE2(u_stage1, u_cur, rhs, dt)
%
% Arguments:
% - u_stage1:    u_approx(t_cur+dt)
% - u_cur:       u(t_cur)
% - rhs:         right-hand side of time evolution equation at (t_cur+dt)
% - dt:          step size
%
% Return value:
% - u_stage2:    u(t_cur + dt)
%
% NOTES:
% - u_stage1, u_cur and rhs are assumed to have the same dimensions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright:  (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:   $Revision$
% Modified:   $Date$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_stage2 = TVD_RK3_STAGE2(u_stage1, u_cur, rhs, dt)

u_stage2 = 0.75*u_cur + 0.25*(u_stage1 + dt*rhs);
