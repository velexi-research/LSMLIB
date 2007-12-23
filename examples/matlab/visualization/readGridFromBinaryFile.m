%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function grid = readGridFromBinaryFile(file_name)
%
% Arguments
% file_name - file name 
% 
% Returns
% grid      - Grid structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright:  (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
% Revision:   $Revision: 1.1 $
% Modified:   $Date: 2006/07/07 12:11:12 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function grid = readGridFromBinaryFile(file_name)

fid = fopen(file_name,'r');

grid.num_dims = fread(fid,1,'int');
grid.x_lo = fread(fid,3,'double');
grid.x_hi = fread(fid,3,'double');
grid.x_lo_ghostbox = fread(fid,3,'double');
grid.x_hi_ghostbox = fread(fid,3,'double');
grid.grid_dims = fread(fid,3,'int');
grid.grid_dims_ghostbox = fread(fid,3,'int');
grid.dx = fread(fid,3,'double');
grid.num_gridpts = fread(fid,1,'int');

fclose(fid);
