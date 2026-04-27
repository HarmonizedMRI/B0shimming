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

% create grid
[X,Y,Z] = ndgrid(1:nx, 1:ny, 1:nz);

% center and scale to cm
X = FOV(1)/nx*(X-nx/2-1/2);
Y = FOV(2)/ny*(Y-ny/2-1/2);
Z = FOV(3)/nz*(Z-nz/2-1/2);

