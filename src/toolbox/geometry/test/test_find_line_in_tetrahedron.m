%
% This script is intended to be used in conjunction with the 
% test_find_line_in_tetrahedron* test programs to verify the correctness 
% of the fineLineInTetrahedron function().  
% 
% This script computes the coefficients for the linear approximation for
% phi and psi within a tetrahedron.  It also computes the line corresponding
% to the intersection of the two planes, and checks whether the specified
% points lie on the line.
%
% Usage: 
% - set coordinates of corners of tetrahedron and 
% - set values of phi and psi at corners
% - set points to check
% - run script
%
% Kevin T. Chu
% MAE, Princeton University
% May 2005
%

% clear workspace
clear;
format long;

%%%%%%%%%%%%%%% USER INPUT SECTION %%%%%%%%%%%%%%%%%%%
X(1,1) = -3.0;
X(1,2) = -6.0;
X(1,3) = 2.0;
X(2,1) = 3.5;
X(2,2) = 4.0;
X(2,3) = -2.0;
X(3,1) = 1.5;
X(3,2) = 4.0;
X(3,3) = -2.0;
X(4,1) = -0.75;
X(4,2) = -1.0;
X(4,3) = -2.0;

phi = [1.0; 1.0; -1.0; -2.0];
psi = [-2.0; -2.0; 2.0; -1.0];

points_to_check = [-0.750000,-1.000000,0.000000];

%%%%%%%%%%%%%%%%% MAIN COMPUTATION %%%%%%%%%%%%%%%%%%%

X = [ones(4,1), X];

% compute coefficients of linear approximation for phi and psi
alpha = X\phi;  alpha = alpha'
beta  = X\psi;  beta = beta'

% compute direction of intersection
line_direction = cross(alpha(2:end),beta(2:end));
line_direction = line_direction/line_direction(3)

% find point on line
if ( det([alpha(2:3);beta(2:3)]) ~= 0 )
  xy = [alpha(2:3);beta(2:3)] \ [-alpha(1); -beta(1)];
  point_on_line = [xy' 0]
else
end

% check if points_to_check lie on the line
num_points = size(points_to_check,1);
point_check_results = points_to_check - repmat(point_on_line,num_points,1);
point_check_results = cross(point_check_results, ...
                            repmat(line_direction,num_points,1))
