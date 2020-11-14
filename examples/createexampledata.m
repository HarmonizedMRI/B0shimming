% Generate example data for testing B0 shimming optimization tools
%
% Output:
%   example.mat, containing:
%   A     9x9            calibration matrix (scanner-specific units)
%   f0    [nx ny nz]     noisy fieldmap (Hz)
%   X     [nx ny nz]     Voxel x coordinate (cm)
%   Y     [nx ny nz]     Voxel y coordinate (cm)
%   Z     [nx ny nz]     Voxel z coordinate (cm)
%   mask  [nx ny nz]     spherical mask (bool)
%   

addpath ..   % path to +shim package

% spherical harmonic basis degree
l = 4;

%% Calculate calibration matrix A
% This only needs to be done once per scanner

% load calibration data
load Ffull                             % [nx ny nz 8]  (does not include DC offset term)
[nx, ny, nz, nShim] = size(Ffull);
fov = 19.8*[1 1 1];   % cm

% generate mask (only for purpose of calculating A)
mask = sum(abs(Ffull),4) > 10;         % [nx ny nz]     
N = sum(mask(:));

% create F = [N 8] by applying mask 
F = zeros(N, nShim);
for ii = 1:nShim
	tmp = Ffull(:,:,:,ii);
	F(:,ii) = tmp(mask); 
end

% specify shim amplitudes (differences) used to obtain F (hardware units)
S = diag([20*ones(1,3) 100*ones(1,5)]);  % Do not include DC term. Hardware units (2nd order unit = 10 mA)

% Get spherical harmonic basis evaluated at the same N spatial locations as F
[X,Y,Z] = shim.getgrid(nx,ny,nz,fov);
H = shim.getSHbasis(X(mask),Y(mask),Z(mask),l);   % [N 2l+1]

% Calculate calibration matrix A
% F and S do NOT include the DC term, i.e., F: [N 8], S: [8 8]
A = shim.getcalmatrix(F, H, S);


%% synthesize noisy fieldmap f0 
nx = 60; ny = 60; nz = 20;
fov = [20 20 20];   % cm

[X,Y,Z] = shim.getgrid(nx,ny,nz,fov);
R = sqrt(X.^2 + Y.^2 + Z.^2);
mask = ones(size(X));
mask(R > 0.8*fov(1)/2) = 0;
mask = logical(mask);
H = shim.getSHbasis(X(mask),Y(mask),Z(mask),l);   % [N 2l+1]
s(1:4) = 10*randn([1 4]);
s(5:9) = 20*randn([1 5]);
R = sqrt(abs(X(mask).^2+Y(mask).^1.5)); 
f0 = 10*H*A*s(:)./(10+R.^1.5);
f0 = f0 + randn(size(f0))*max(f0(:))/10;
f0 = embed(f0, mask);

%% Save to file
save exampledata A f0 X Y Z mask



