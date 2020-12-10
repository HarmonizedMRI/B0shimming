function dat = loadScanArchive(fname)
% function dat = loadScanArchive(fname)
%
% Output:
%  dat  [ndat ncoils nframes]

% addpath /mnt/storage/jfnielse/Orchestra/orchestra-sdk-1.9-1.matlab/   

archive = GERecon('Archive.Load', fname);   %'ScanArchive_7347633TMRFIX_20201010_170753212.h5');

if nargin < 2
	frameStart = 1;
	frameStop = archive.FrameCount;
end

% initialize data output matrix
currentControl = GERecon('Archive.Next', archive);
tmp = currentControl.Data;   % size = [ndat ncoils]
[ndat ncoils] = size(tmp);
fprintf('Initializing output data matrix, size [%d %d %d]\n', ndat, ncoils, archive.FrameCount);
dat = zeros(ndat, ncoils, archive.FrameCount);

% load data
dat(:,:,1) = currentControl.Data; 
for ii = 2:archive.FrameCount
	if ~mod(ii,1000)
		for ibs = 1:30; fprintf('\b'); end;
		fprintf('Frame %d of %d', ii, archive.FrameCount);
	end
	currentControl = GERecon('Archive.Next', archive);
	dat(:,:,ii) = currentControl.Data;   % size = [ndat ncoils nframes]
end

return

