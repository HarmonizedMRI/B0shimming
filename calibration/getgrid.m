function [X,Y,Z] = getgrid(nx,ny,nz,FOV)
% function [X,Y,Z] = getgrid(nx,ny,nz,FOV)
%
% Create 3D grid in distance units (cm), centered about 0.
%
% Inputs
%  nx, ny, nz    int       grid size
%  FOV           [3 1]     field of view in x/y/z (cm)
%
% Outputs
%  X     [nx ny nz]    x coordinate (cm)
%  Y     [nx ny nz]    y coordinate
%  Z     [nx ny nz]    z coordinate

% normalized grid
[X,Y,Z] = ndgrid(linspace(-1+2/nx/2, 1-2/nx/2, nx), ...
   linspace(-1+2/ny/2, 1-2/ny/2, ny), ...
	linspace(-1+2/nz/2, 1-2/nz/2 ,nz));

% scale to cm
X = FOV(1)*X;
Y = FOV(2)*Y;
Z = FOV(3)*Z;

return
