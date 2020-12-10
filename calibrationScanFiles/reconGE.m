function reconGE(dat, readoutFile, imsize, deltaTE)
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

for ie = 1:length(deltaTE)
	offset = (ie-1)*ny*nz;
	d = dat(:, :, (offset+1):(offset+nz*ny) );
	d = reshape(d, [ndat ncoils ny nz]);
	d = permute(d, [1 3 4 2]);                 % [ndat ny nz ncoils]

	for c = 1:ncoils
		imtmp = fftshift(ifftn(fftshift(d(:,:,:,c))));
		ims(:,:,:,c) = imtmp((2*nx-nx/2):(2*nx+nx/2-1), :, :);
	end
	imsos(:,:,:,ie) = sqrt(sum(abs(ims).^2,4));
end

