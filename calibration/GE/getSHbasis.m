function H = getSHbasis(X,Y,Z)
% function H = getSHbasis(X,Y,Z)
%
% get spherical harmonics basis
%
% Inputs:
%  X/Y/Z    [nx ny nz]    x/y/z coordinate corresponding to 3D image matrix
%
% Output:
%   H       [nx*ny*nz 9]    SH basis, including dc (B0 offset) term

x = X(:); y = Y(:); z = Z(:);

% from shimcal.pl: shims = ('x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy');  
H = [ones(numel(X),1) x y z z.^2 x.*y z.*x x.^2-y.^2 z.*y];

return;
