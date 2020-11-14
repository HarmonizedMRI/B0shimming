% demoWLS.m
%
% Demonstrate the proposed B0 shimming workflow using a synthesized baseline fieldmap f0.
% 
% Goal is to calculate shim setting changes s that minimize the weighted RMS B0 inhomogeneity.
% Provided here as an example; other loss functions may be more useful depending on the application.

% spherical harmonic basis degree
l = 6;

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
H = shim.getSHbasis(X(mask),Y(mask),Z(mask),l);   % [N ...]

% Calculate calibration matrix A
% F and S do NOT include the DC term, i.e., F: [N 8], S: [8 8]
A = shim.getcalmatrix(F, H, S);


%% get shim amplitudes 'shat' that minimize RMS fieldmap
%shimNum = 5;
f0 = Ffull(:,:,:,shimNum);  % see how well one shim field fits basis
mask = sum(abs(Ffull),4) > 10;         % [nx ny nz]
f0 = f0(mask);   % [N 1]
W = diag_sp(ones(N,1));
shat = -(W*H*A)\(W*f0);    % [9 1]. NB! May need to be rounded before applying settings on scanner.

f = f0 + H*A*shat;  % predicted fieldmap after applying shims

f0 = embed(f0, mask);
f = embed(f, mask);

% compare baseline and predicted fieldmaps
subplot(121)
im(f0); 
%title(sprintf('fieldmaps: original (left in each panel),\nand after 2nd order shimming (right)'));
title(sprintf('before shimming'));
h = colorbar; h.TickLabels{end} = 'Hz'; % h.Label.String = 'B0 field (Hz)';
subplot(122)
im(f,[-5 5]); 
title(sprintf('Residual (shim term %d)', shimNum));
h = colorbar; h.TickLabels{end} = 'Hz';


