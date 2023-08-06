fn.baseline.a = 'ScanArchive_7347633TMRFIX_20201210_163100493.h5';
fn.z2.plus = 'ScanArchive_7347633TMRFIX_20201210_163254499.h5';
fn.z2.minus = 'ScanArchive_7347633TMRFIX_20201210_163449189.h5';
fn.xy.plus = 'ScanArchive_7347633TMRFIX_20201210_163659066.h5';
fn.xy.minus = 'ScanArchive_7347633TMRFIX_20201210_164001060.h5';
fn.zx.plus = 'ScanArchive_7347633TMRFIX_20201210_164229588.h5';
fn.zx.minus = 'ScanArchive_7347633TMRFIX_20201210_164506004.h5';
fn.x2y2.plus = 'ScanArchive_7347633TMRFIX_20201210_165556534.h5';
fn.x2y2.minus = 'ScanArchive_7347633TMRFIX_20201210_165745213.h5';
fn.zy.plus = 'ScanArchive_7347633TMRFIX_20201210_170051019.h5';
fn.zy.minus = 'ScanArchive_7347633TMRFIX_20201210_170423495.h5';
fn.x.plus = 'ScanArchive_7347633TMRFIX_20201210_170631453.h5';
fn.x.minus = 'ScanArchive_7347633TMRFIX_20201210_171520671.h5';
fn.y.plus = 'ScanArchive_7347633TMRFIX_20201210_171713560.h5';
fn.y.minus = 'ScanArchive_7347633TMRFIX_20201210_171900191.h5';
fn.z.plus = 'ScanArchive_7347633TMRFIX_20201210_172059848.h5';
fn.z.minus = 'ScanArchive_7347633TMRFIX_20201210_172248592.h5';
fn.baseline.b = 'ScanArchive_7347633TMRFIX_20201210_172449859.h5';

shims = {'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'};

% get mask
dat = loadScanArchive(fn.baseline.a);
[img, mask] = reconGE(dat, 'readout.mod', [60 60 50], [0 0.5 1 2]);
%pc = phasecontrastmulticoil(img(:,:,:,:,1), img(:,:,:,:,3));
%pc(~mask) = 0;
dte = 1e-3;   % sec
%f0 = pc/(2*pi*dte);  % Hz

for ii = 1:length(shims)
	dat = loadScanArchive(eval(sprintf('fn.%s.plus', shims{ii})));
	img = reconGE(dat, 'readout.mod', [60 60 50], [0 0.5 1 2]);
	pc = phasecontrastmulticoil(img(:,:,:,:,1), img(:,:,:,:,3));
	fplus = pc/(2*pi*dte);  % Hz
	dat = loadScanArchive(eval(sprintf('fn.%s.minus', shims{ii})));
	img = reconGE(dat, 'readout.mod', [60 60 50], [0 0.5 1 2]);
	pc = phasecontrastmulticoil(img(:,:,:,:,1), img(:,:,:,:,3));
	fminus = pc/(2*pi*dte);  % Hz
	f = fplus-fminus;
	f(~mask) = 0;
	F(:,:,:,ii) = f;
	subplot(3,3,ii); im(F(:,:,:,ii)); title(shims{ii}); colormap jet;
	pause(1);
end	

return;

dat = loadScanArchive(fn.xy.plus);
img = reconGE(dat, 'readout.mod', [60 60 50], [0 0.5 1 2]);
pc = phasecontrastmulticoil(img(:,:,:,:,1), img(:,:,:,:,3));
pc(~mask) = 0;
fxy = pc/(2*pi*dte);  % Hz
figure; im(fxy-f0); title('xy plus'); colormap jet;  % l,m = 2,-2

dat = loadScanArchive(fn.zy.plus);
img = reconGE(dat, 'readout.mod', [60 60 50], [0 0.5 1 2]);
pc = phasecontrastmulticoil(img(:,:,:,:,1), img(:,:,:,:,3));
pc(~mask) = 0;
fzy = pc/(2*pi*dte);  % Hz
figure; im(fzy-f0); title('zy plus'); colormap jet;  % l,m = 2,-1 

dat = loadScanArchive(fn.x2y2.plus);
img = reconGE(dat, 'readout.mod', [60 60 50], [0 0.5 1 2]);
pc = phasecontrastmulticoil(img(:,:,:,:,1), img(:,:,:,:,3));
pc(~mask) = 0;
fx2y2 = pc/(2*pi*dte);  % Hz
figure; im(fx2y2-f0); title('x2y2 plus'); colormap jet;  % l,m = 2,2


nx = size(img,1); ny = size(img,2); nz = size(img,3);
res = 24/60;   % voxel size (isotropic)
xrange = res*[(-nx/2+1/2):(nx/2-1/2)];  % may be 1/2 pixel off here...
yrange = res*[(-ny/2+1/2):(ny/2-1/2)];
zrange = res*[(-nz/2+1/2):(nz/2-1/2)];
[x,y,z] = ndgrid(xrange, yrange, zrange);

% display spherical harmonics for comparison
if false
for l = 0:2
	for m = -l:l
		fsh = evalspharm(x(:), y(:), z(:), l, m);
		fsh = reshape(fsh, [nx ny nz]);
		figure; im(real(fsh)); title(sprintf('l, m = %d, %d', l, m)); colormap jet;
		figure; im(imag(fsh)); title(sprintf('l, m = %d, %d', l, -m)); colormap jet;
		
	end
end
end


