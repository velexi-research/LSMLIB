%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% computePlusAndMinusDerivatives1D() computes the upwind Hamilton-Jacobi
% spatial derivatives of the specified order for 1D grid functions.
%
% Usage: [phi_x_plus, phi_x_minus] = ...
%        computePlusAndMinusDerivatives1D(phi, ...
%                                         ghostcell_width, ...
%                                         dx, ...
%                                         spatial_derivative_order)
%
% Arguments:
% - phi:                       function for which to compute upwind 
%                                derivative
% - ghostcell_width:           number of ghostcells at boundary of 
%                                computational domain
% - dx:                        grid cell size
% - spatial_derivative_order:  spatial derivative order
%                                (default = 1)
%
% Return values:
% - phi_x_plus:                plus HJ ENO/WENO derivative
% - phi_x_minus:               minus HJ ENO/WENO derivative
%
% NOTES:
% - phi_x_plus and phi_x_minus have the same ghostcell width as phi.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright:  (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
% Revision:   $Revision: 1.1 $
% Modified:   $Date$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phi_x_plus, phi_x_minus] = ...
         computePlusAndMinusDerivatives1D(phi, ...
                                          ghostcell_width, ...
                                          dx, ...
                                          spatial_derivative_order)
 

% parameter checks
if (nargin < 3)
  error('MATLAB:missingArgs','computePlusAndMinusDerivatives1D:missing arguments');
end

if (nargin < 4)
  spatial_derivative_order = 1;
else
  if ( (spatial_derivative_order ~= 1) & (spatial_derivative_order ~= 2) ...
     & (spatial_derivative_order ~= 3) & (spatial_derivative_order ~= 5) )

    error('computePlusAndMinusDerivatives1D:Invalid spatial derivative order...only 1, 2, 3, and 5 are supported');
  end
end

switch (spatial_derivative_order)

  case 1
    [phi_x_plus, phi_x_minus] = HJ_ENO1_1D(phi, ghostcell_width, dx);

  case 2
    [phi_x_plus, phi_x_minus] = HJ_ENO2_1D(phi, ghostcell_width, dx);

  case 3
    [phi_x_plus, phi_x_minus] = HJ_ENO3_1D(phi, ghostcell_width, dx);

  case 5
    [phi_x_plus, phi_x_minus] = HJ_WENO5_1D(phi, ghostcell_width, dx);

end
