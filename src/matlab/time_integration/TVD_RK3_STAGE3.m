%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TVD_RK3_STAGE3 advances the solution through the final stage of a 
%                third-order TVD Runge-Kutta step.
%
% Usage: u_next = TVD_RK3_STAGE3(u_stage2, u_cur, rhs, dt)
%
% Arguments:
% - u_stage2:    u_approx(t_cur+dt)
% - u_cur:       u(t_cur)
% - rhs:         right-hand side of time evolution equation at (t_cur+dt/2)
% - dt:          step size
%
% Return value:
% - u_next:      u(t_cur + dt)
%
% NOTES:
% - u_stage2, u_cur and rhs are assumed to have the same dimensions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright:  (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:   $Revision: 1.3 $
% Modified:   $Date: 2006/04/22 12:39:34 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_next = TVD_RK3_STAGE3(u_stage2, u_cur, rhs, dt)

u_next = u_cur/3.0 + 2/3*(u_stage2 + dt*rhs);
