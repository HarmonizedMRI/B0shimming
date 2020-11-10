% Generate example data for testing B0 shimming optimization tools
%
% Output:
%   example.mat, containing:
%   A   9x9     calibration matrix
%   f0  [N 1]   noisy fieldmap
%   X
%   

addpath ..   % path to +shim package


%% Calculate calibration matrix A
% This only needs to be done once per scanner

% load calibration data
load Ffull                             % [nx ny nz 8]  (does not include DC offset term)
[nx, ny, nz, nShim] = size(Ffull);
FOV = 19.8*[1 1 1];   % cm

% generate mask
mask = sum(abs(Ffull),4) > 10;         % [nx ny nz]     
N = sum(mask(:));

% create F by applying mask 
F = zeros(N, nShim);
for ii = 1:nShim
	tmp = Ffull(:,:,:,ii);
	F(:,ii) = tmp(mask); 
end

% specify shim amplitudes (differences) used to obtain F (hardware units)
S = diag([20*ones(1,3) 100*ones(1,5)]);  % do not include DC term. 2nd order amp units = 10 mA

% Get spherical harmonic basis evaluated at the same N spatial locations as F
[X,Y,Z] = shim.getgrid(nx,ny,nz,FOV);
H = shim.getSHbasis(X(mask),Y(mask),Z(mask));   % [N 9]

% Calculate calibration matrix A
% F and S do NOT include the DC term, i.e., F: [N 8], S: [8 8]
A = shim.getcalmatrix(F, H, S)

%% synthesize noisy fieldmap f0 = [N 1]
nx = 60; ny = 60; nz = 20;
FOV = [20 20 8];   % cm

[X,Y,Z] = shim.getgrid(nx,ny,nz,FOV);
R = sqrt(X.^2 + Y.^2 + Z.^2);
m = ones(size(X));
m(R < 0.8*FOV/2
H = shim.getSHbasis(X(:),Y(:),Z(:));   % [N 9]
s(1:4) = 10*randn([1 4]);
s(5:9) = 20*randn([1 5]);
R = sqrt(X(:).^2+Y(:).^1.5);
f0 = 10*H*A*s(:)./(10+R.^1.5);
f0 = f0 + randn(size(f0))*max(f0(:))/10;

%% Save A and f0 to file
% get shim amplitudes 'shat' that minimize RMS fieldmap
N = size(f0,1);
W = diag_sp(ones(N,1));
shat = -(W*H*A)\(W*f0(:));    % [9 1]. NB! May need to be rounded before applying settings on scanner.

% compare baseline and predicted fieldmaps
f = f0 + H*A*shat;
f = reshape(f, [nx ny nz]);
f0 = reshape(f0, [nx ny nz]);
im(cat(1, f0, f)); 
title(sprintf('fieldmaps: original (left in each panel),\nand after 2nd order shimming (right)'));
h = colorbar; h.TickLabels{end} = 'Hz'; % h.Label.String = 'B0 field (Hz)';



