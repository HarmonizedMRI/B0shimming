b0 = reconB0(d, deltaTE)
% b0 = reconB0(d, deltaTE)
%
% Reconstruct B0 field map from date obtained with b0.seq, 
% see writeB0.m.
%
% Inputs:
%   dat          [nx ny nz ncoils 2]    Complex multicoil data at two echo times
%   deltaTE      [1]                    TE difference (sec)                    
%
% Outputs:
%   b0           [nx ny nz]             field map (Hz)
%
% Example:
%   >> deltaTE = 2.2369e-3;  % TE difference at 3.0T. See writeB0.m
%   >> [b0, mask] = reconB0(dat, deltaTE);
%   >> im(b0.*mask); colormap parula;  % Requires MIRT toolbox, https://github.com/JeffFessler/mirt
%   >> colorbar; h.TickLabels(end) = {'Hz'};

% reconstruct complex coil images
d_te1 = d(:,:,:,1);   % TE1 data, size [60 60 60]
d_te2 = d(:,:,:,2);   % TE1 data, size [60 60 60]
for ic = 1:ncoil
    ims_te1(:,:,:,ic) = fftshift(ifftn(fftshift(d_te1(:,:,:,ic))));
    ims_te2(:,:,:,ic) = fftshift(ifftn(fftshift(d_te2(:,:,:,ic))));
end

% phase difference map, radians
pc = toppe.utils.phasecontrastmulticoil(ims_te2, ims_te1);

% field map in Hz
b0 = pc/(2*pi)/deltaTE;

% object support
I_te1 = sqrt(sum(abs(ims_te1).^2, ndims(d)));   % root sum of squares coil combination
mask = I_te1 > 0.1*max(I_te1(:));
