%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% computeUpwindDerivatives1D() computes the upwind Hamilton-Jacobi
% spatial derivatives of the specified order for 1D grid functions.
%
% Usage: phi_x = computeUpwindDerivatives1D(phi, ...
%                                           vel_x, ...
%                                           ghostcell_width, ...
%                                           dx, ...
%                                           spatial_derivative_order)
%
% Arguments:
% - phi:                       function for which to compute upwind 
%                                derivative
% - vel_x:                     velocity to use in upwinding
% - ghostcell_width:           number of ghostcells at boundary of 
%                                computational domain
% - dx:                        grid cell size
% - spatial_derivative_order:  spatial derivative order
%                                (default = 1)
%
% Return values:
% - phi_x:                     upwind HJ ENO/WENO derivative
%
% NOTES:
% - phi_x has the same ghostcell width as phi.
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

function phi_x = computeUpwindDerivatives1D(phi, ...
                                            vel_x, ...
                                            ghostcell_width, ...
                                            dx, ...
                                            spatial_derivative_order)

% parameter checks
if (nargin < 4)
  error('MATLAB:missingArgs','computeUpwindDerivatives1D:missing arguments');
end

if (nargin < 5)
  spatial_derivative_order = 1;
else
  if ( (spatial_derivative_order ~= 1) & (spatial_derivative_order ~= 2) ...
     & (spatial_derivative_order ~= 3) & (spatial_derivative_order ~= 5) )

    error('computeUpwindDerivatives1D:Invalid spatial derivative order...only 1, 2, 3, and 5 are supported');
  end
end

switch (spatial_derivative_order)

  case 1
    phi_x = UPWIND_HJ_ENO1_1D(phi, vel_x, ghostcell_width, dx);

  case 2
    phi_x = UPWIND_HJ_ENO2_1D(phi, vel_x, ghostcell_width, dx);

  case 3
    phi_x = UPWIND_HJ_ENO3_1D(phi, vel_x, ghostcell_width, dx);

  case 5
    phi_x = UPWIND_HJ_WENO5_1D(phi, vel_x, ghostcell_width, dx);

end
