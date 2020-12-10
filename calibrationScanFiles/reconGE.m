function F = reconGE(dat, readoutFile, imSize, deltaTE)
%
% Example usage:
% >> dat = loadScanArchive('ScanArchive_7347633TMRFIX_20201210_120317865.h5');
% >> readoutFile = 'readout.mod';
% >> imSize = [60 60 50];
% >> deltaTE = [0 0.5 1.0 2.0];    % ms


nx = imSize(1);
ny = imSize(2);
nz = imSize(3);

% crop out readout plateau
[~,~,~,~,~,paramsint16] = toppe.readmod(readoutFile);
daqStart = paramsint16(1) + 1;
daqStop = daqStart + paramsint16(2)-1;
dat = dat(daqStart:daqStop, :, :, :);

ndat = size(dat, 1);
ncoils = size(dat, 2);
nframes = size(dat, 3);
nechoes = length(deltaTE);

d = zeros(ndat, ncoils, ny, nz, nechoes);

ifr = 1;    % frame number
for iz = 1:nz
	for iy = 1:ny
		for ie = 1:nechoes
			d(:,:,iy,iz,ie) = dat(:,:,ifr);
			ifr = ifr + 1;
		end
	end
end

d = permute(d, [1 3 4 2 5]);      % [ndat ny nz ncoils nechoes]

ims = zeros(nx, ny, nz, ncoils, nechoes);
for ie = 1:nechoes
	for c = 1:ncoils
		imtmp = fftshift(ifftn(fftshift(d(:,:,:,c,ie))));
		ims(:,:,:,c,ie) = imtmp((2*nx-nx/2):(2*nx+nx/2-1), :, :);
	end
end

F = ims;
