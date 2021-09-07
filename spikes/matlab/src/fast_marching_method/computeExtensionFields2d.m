%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% computeExtensionFields2d() computes a distance from an 
% arbitrary level set function using the Fast Marching Method.
% 
% Usage: [distance_function, extension_fields] = ...
%        computeExtensionFields2d(phi, source_fields, ...
%                                 dX, mask, spatial_discretization_order)
%
% Arguments:
% - phi:                           level set function to use in 
%                                    computing distance function
% - source_fields:                 field variables that are to
%                                    be extended off of the zero
%                                    level set
% - dX:                            array containing the grid spacing
%                                    in each coordinate direction
% - mask:                          mask for domain of problem;
%                                    grid points outside of the domain
%                                    of the problem should be set to a
%                                    negative value
%                                    (default = [])
% - spatial_discretization_order:  order of discretization for 
%                                    spatial derivatives
%                                    (default = 2)
%
% Return values:
% - distance_function:             distance function
% - extension_fields:              extension fields
%
% NOTES:
% - All data arrays are assumed to be in the order generated by the 
%   MATLAB meshgrid() function.  That is, data corresponding to the 
%   point (x_i,y_j) is stored at index (j,i) in the data array.
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
