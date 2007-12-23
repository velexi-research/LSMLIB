%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [data,nx,ny,nz] = readDataArray(file_name)
%
% Reads a binary data file which contains 3 binary integers nx, ny, nz 
% indicating sizes of the array in each dimension, followed by an array
% of nx*ny*nz double precision real numbers.
% 
% Arguments:
% file_name - file name 
%
% Returns
% data      - double array read from the file
% nx,ny,nz  - array dimensions 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright:  (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
% Revision:   $Revision: 1.1 $
% Modified:   $Date: 2006/07/07 12:11:12 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,nx,ny,nz] = readDataArray(file_name)

fid = fopen(file_name,'r');

nx = fread(fid,1,'int');
ny = fread(fid,1,'int');
nz = fread(fid,1,'int');

if(nz > 1)
    data = fread(fid,nx*ny*nz,'double');
    data = reshape(data,[nx ny nz]);
else
    data = fread(fid,nx*ny,'double');
    data = reshape(data,[nx ny]);
end

fclose(fid);
