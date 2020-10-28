function [s,f] = getshimsWLS(f0, H, A, W)
% function [s,f] = getshimsWLS(f0, H, A, W = diag_sp(ones(N,1)))
%
% Return shim setting changes s that minimize the weighted RMS B0 inhomogeneity.
% Provided here as an example; other loss functions may be more useful depending on the application.
%
% Inputs:
%   f0     [N 1]      (Hz) Phase-unwrapped 3D fieldmap
%   H      [N 9]      SH basis matrix, evaluated at the same spatial locations (control points) as f0.
%   A      [9 9]      shim calibration matrix. See getcalmatrix.m
%   W      [N N]      optional sparse diagonal weighting matrix. Default = diag_sp(ones(N,1)) (from MIRT toolbox)
% 
% Output:
%   s      [9 1]      change in shim settings from the baseline settings used to acquire f0
%                     s(1) is the B0 frequency offset (Hz)
%                     s(2:9) are the linear and 2nd order shim amplitude changes in hardware units (e.g., mA)
%   f      [N 1]      predicted fieldmap after applying s (f = H*A*s + f0)

N = size(f0,1);

if nargin < 4
	W = diag_sp(ones(N,1));
end

s = -(W*H*A)\(W*f0(:));    % [9 1]. NB! May need to be rounded before applying settings on scanner.

f = f0 + H*A*s;

return;

% display predicted f 
% display transpose of images, i.e., frequency-encoding direction is horizontal, to match GUI
iz = izDisplay; %nz/2-10;
f2d = fmap(:,:,iz);
fmin = min(f2d(:));
fmax = max(f2d(:));

fmin = -100;
fmax = 100;

subplot(131); imagesc(fmap(:,:,iz)'); title('acquired field map'); 
axis equal off; caxis([fmin fmax]/2);
xlabel('L                               R');
ylabel('P                               A');
h = colorbar; title(h, 'Hz');

fmap_applied = reshape((C*A)*s, [nx ny nz]).*mask; 
%fmap_applied = fmap_applied - mean(fmap_applied(logical(mask)));

subplot(132); imagesc(-fmap_applied(:,:,iz)'); title('spherical harmonic fit'); 
axis equal off; caxis([fmin fmax]/2);
h = colorbar; title(h, 'Hz');

fmap = fmap.*mask;
subplot(133); imagesc(fmap(:,:,iz)'+fmap_applied(:,:,iz)'); title('after applying new shims');
axis equal off; caxis([fmin fmax]/2);
h = colorbar; title(h, 'Hz');

return
