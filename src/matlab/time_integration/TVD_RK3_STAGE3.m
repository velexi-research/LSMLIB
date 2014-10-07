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
% Copyrights: (c) 2005 The Trustees of Princeton University and Board of
%                 Regents of the University of Texas.  All rights reserved.
%             (c) 2009 Kevin T. Chu.  All rights reserved.
% Revision:   $Revision$
% Modified:   $Date$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_next = TVD_RK3_STAGE3(u_stage2, u_cur, rhs, dt)

u_next = u_cur/3.0 + 2/3*(u_stage2 + dt*rhs);
