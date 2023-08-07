b0 = reconB0(d, deltaTE, thresh)
% b0 = reconB0(d, deltaTE, [thresh=0.1])
%
% Reconstruct B0 field map from date obtained with b0.seq, 
% see writeB0.m.
%
% Inputs:
%   d            [nx ny nz ncoils 2]    Complex multicoil data at two echo times
%   deltaTE      [1]                    TE difference (sec)                    
%   thresh       [1]                    Magnitude threshold for defining mask. Default = 0.1
%
% Outputs:
%   b0           [nx ny nz]             field map (Hz)
%
% Example:
%   >> dat = loaddata_ge('P12345.7');
%   >> deltaTE = 2.2369e-3;  % TE difference at 3.0T. See writeB0.m
%   >> [b0, mask] = reconB0(dat, deltaTE);
%   >> im(b0.*mask); colormap parula;  % Requires MIRT toolbox, https://github.com/JeffFessler/mirt
%   >> colorbar; h.TickLabels(end) = {'Hz'};

if nargin < 3
    thresh = 0.1;
end

% reconstruct complex coil images
d_te1 = d(:,:,:,:,1);   % TE1 data
d_te2 = d(:,:,:,:,2);   % TE2 data
for ic = 1:ncoil
    ims_te1(:,:,:,ic) = fftshift(ifftn(fftshift(d_te1(:,:,:,ic))));
    ims_te2(:,:,:,ic) = fftshift(ifftn(fftshift(d_te2(:,:,:,ic))));
end

% phase difference map, radians. Magnitude-squared coil weighting.
pc = toppe.utils.phasecontrastmulticoil(ims_te2, ims_te1);

% field map in Hz
b0 = pc/(2*pi)/deltaTE;

% object support
I_te1 = sqrt(sum(abs(ims_te1).^2, ndims(d)));   % root sum of squares coil combination
mask = I_te1 > thresh*max(I_te1(:));
