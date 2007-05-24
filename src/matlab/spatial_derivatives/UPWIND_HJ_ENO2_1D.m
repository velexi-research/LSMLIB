%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% UPWIND_HJ_ENO2_1D() computes the second-order upwind HJ ENO
% approximation to phi_x.
%
% Usage: phi_x = UPWIND_HJ_ENO2_1D(phi, vel_x, ghostcell_width, dx)
%
% Arguments:
% - phi:               function for which to compute upwind 
%                        derivative
% - vel_x:             velocity to use in upwinding
% - ghostcell_width:   number of ghostcells at boundary of 
%                        computational domain
% - dx:                grid cell size
%
% Return values:
% - phi_x:             second-order, upwind HJ ENO derivative
%
% NOTES:
% - phi_x has the same ghostcell width as phi.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author:     Kevin T. Chu
% Copyright:  (c) 2005-2006, MAE Princeton University 
% Revision:   $Revision: 1.3 $
% Modified:   $Date: 2006/04/24 00:52:39 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
