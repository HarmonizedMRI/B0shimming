% demoWLS.m
%
% Demonstrate the proposed B0 shimming workflow using a synthesized baseline fieldmap f0.
% 
% Goal is to calculate shim setting changes s that minimize the weighted RMS B0 inhomogeneity.
% Provided here as an example; other loss functions may be more useful depending on the application.
%
% Inputs:
%   f0     [N 1]      (Hz) Phase-unwrapped 3D fieldmap. N = number of voxels ('control points').
%   H      [N 9]      SH basis matrix, evaluated at the same spatial locations (control points) as f0.
%   A      [9 9]      shim calibration matrix. See getcalmatrix.m
%   W      [N N]      optional sparse diagonal weighting matrix. Default = diag_sp(ones(N,1)) (from MIRT toolbox)
% 
% Output:
%   s      [9 1]      change in shim settings from the baseline settings used to acquire f0
%                     s(1) is the B0 frequency offset (Hz)
%                     s(2:9) are the linear and 2nd order shim amplitude changes in hardware units (e.g., mA)
%   f      [N 1]      predicted fieldmap after applying s (f = H*A*s + f0)

addpath ..   % path to +shim package

% load calibration matrix
load A

% synthesize noisy fieldmap f0 = [N 1]
nx = 64; ny = 64; nz = 9;
FOV = [20 20 5];
[X,Y,Z] = shim.getgrid(nx,ny,nz,FOV);
H = shim.getSHbasis(X(:),Y(:),Z(:));   % [N 9]
struth(1:4) = 0.2*randn([1 4]);
s(5:9) = 1e0*randn([1 5]);
f0 = reshape(H*A*s(:), [nx ny nz]);
f0 = f0 + randn(size(f0))*max(f0(:))/10;
f0 = f0(:);

% get shim adjustment that minimizes RMS fieldmap
N = size(f0,1);
W = diag_sp(ones(N,1));
shat = -(W*H*A)\(W*f0(:));    % [9 1]. NB! May need to be rounded before applying settings on scanner.

% compare baseline and predicted fieldmaps
f = f0 + H*A*shat;
f = reshape(f, [nx ny nz]);
f0 = reshape(f0, [nx ny nz]);

im(cat(1, f0, f)); colorbar;

return;

