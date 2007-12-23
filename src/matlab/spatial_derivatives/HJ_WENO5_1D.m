%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HJ_WENO5_1D() computes the fifth-order plus and minus HJ WENO
% approximation to phi_x.
%
% Usage: [phi_x_plus, phi_x_minus] = HJ_WENO5_1D(phi, ghostcell_width, dx)
%
% Arguments:
% - phi:               function for which to compute plus and minus
%                        spatial derivatives
% - ghostcell_width:   number of ghostcells at boundary of 
%                        computational domain
% - dx:                grid cell size
%
% Return values:
% - phi_x_plus:        fifth-order, plus HJ WENO derivative
% - phi_x_minus:       fifth-order, minus HJ WENO derivative
%
% NOTES:
% - phi_x_plus and phi_x_minus have the same ghostcell width as phi.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright:  (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:   $Revision: 1.3 $
% Modified:   $Date: 2006/04/24 00:52:38 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
