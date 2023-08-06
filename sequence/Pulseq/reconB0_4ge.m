% reconB0.m
%
% Reconstruct B0 field maps obtained with b0.seq

% load data from P-file (set CV2 = 1 to save P-file)
d = toppe.utils.loadpfile('P,b0.7');
[nfid,ncoil,nslice,necho,nview] = size(d);   % necho = 1 by construction
d = permute(d, [1 5 3 2 4]);   % [nfid nview nslice ncoil]

% discard data during receive gain calibration (see writeB0.m)
d = d(:,:,2:end,:);  

% reconstruct complex coil images
d_te1 = d(:,1:2:end,:,:);   % TE1 data, size [60 60 60]
d_te2 = d(:,2:2:end,:,:);   % TE2 data, size [60 60 60]
for ic = 1:ncoil
    ims_te1(:,:,:,ic) = fftshift(ifftn(fftshift(d_te1(:,:,:,ic))));
    ims_te2(:,:,:,ic) = fftshift(ifftn(fftshift(d_te2(:,:,:,ic))));
end

% phase difference map, radians
pc = toppe.utils.phasecontrastmulticoil(ims_te2, ims_te1);

% TE difference at 3.0 T, see writeB0.m
dte = 2.2369e-3;  % sec

% field map (Hz)
b0 = pc/(2*pi)/dte;

% Display. Requires MIRT toolbox, https://github.com/JeffFessler/mirt
I_te1 = sqrt(sum(abs(ims_te1).^2, ndims(d)));   % root sum of squares coil combination
mask = I_te1 > 0.1*max(I_te1(:));
figure; im(b0.*mask); colormap parula;
h = colorbar; h.TickLabels(end) = {'Hz'};
