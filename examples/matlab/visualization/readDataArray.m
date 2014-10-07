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
% Copyrights: (c) 2005 The Trustees of Princeton University and Board of
%                 Regents of the University of Texas.  All rights reserved.
%             (c) 2009 Kevin T. Chu.  All rights reserved.
% Revision:   $Revision$
% Modified:   $Date$
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
