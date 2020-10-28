function [X,Y,Z] = getgrid(nx,ny,nz)
% get grid coordinates normalized to (-1, 1), symmetric about zero

[X,Y,Z] = ndgrid(linspace(-1+2/nx/2, 1-2/nx/2, nx), ...
   linspace(-1+2/ny/2, 1-2/ny/2, ny), ...
	linspace(-1+2/nz/2, 1-2/nz/2 ,nz));

return
