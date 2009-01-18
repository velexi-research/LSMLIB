%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function plotMultipleZeroLevelSets(data_filename1,color1,..., data_filenameN,colorN,grid_filename)
%
% Arguments:
% data_filename1,...data_filenameN - file names for data arrays
% color1,...,colorN - color codes (standard Matlab) for each data array
% grid_filename - file name for binary Grid array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright:  (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
% Revision:   $Revision$
% Modified:   $Date$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMultipleZeroLevelSets(varargin)

N = nargin;
M = N-1;

grid_filename = varargin{N};
grid = readGridFromBinaryFile(grid_filename);

% read pairs of data filenames and appropriate colors
j=1;
for i = 1 : 2: M
  data_filename = varargin{i};
  color{j} = varargin{i+1};
  [data_tmp, nx, ny, nz] = readDataArray(data_filename);
  data{j} = data_tmp;
  j = j+1;
end

if( nz == 1) % 2D visualization
    
  x = grid.x_lo_ghostbox(1)*ones(1,nx) + grid.dx(1)*(0:nx-1);
  y = grid.x_lo_ghostbox(2)*ones(1,ny) + grid.dx(2)*(0:ny-1);
  [X Y]= meshgrid(x,y);
  
  figure;
  level = [ 0 0 ];
  
  for i=1:M/2
   data_tmp  = (data{i})'; 
   [ garbage, hI ] = contour(X, Y, data_tmp, level,color{i}); xlabel('x'); ylabel('y');
   daspect([1 1 1]);
   hold on;
  end
  hold off
else % 3D visualization
    
  x = grid.x_lo_ghostbox(1)*ones(1,nx) + grid.dx(1)*(0:nx-1);
  y = grid.x_lo_ghostbox(2)*ones(1,ny) + grid.dx(2)*(0:ny-1);
  z = grid.x_lo_ghostbox(3)*ones(1,nz) + grid.dx(3)*(0:nz-1);
  [Y X Z] = meshgrid(y,x,z);
  
  level = 0;
  
  for i=1:M/2
    h = patch(isosurface(Y,X,Z,data{i},level));
    set(h,'FaceColor',color{i},'EdgeColor', 'none');
    hold on    
  end
  
  daspect([1 1 1]); camlight;
  xlabel('x'); ylabel('y'); zlabel('z');
  view(3)  
  hold off
end    
