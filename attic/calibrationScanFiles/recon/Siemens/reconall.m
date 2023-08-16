fn.baseline = 'img_meas_MID00618_FID31337_pulseq_baseline.mat';
fn.z2.plus = 'img_meas_MID00642_FID31361_pulseq_z2plus20.mat';
fn.z2.minus = 'img_meas_MID00639_FID31358_pulseq_z2minus20.mat';
fn.xy.plus = 'img_meas_MID00660_FID31379_pulseq_xyplus20.mat';
fn.xy.minus = 'img_meas_MID00662_FID31381_pulseq_xyminus20.mat';
fn.zx.plus = 'img_meas_MID00650_FID31369_pulseq_zxplus20.mat';
fn.zx.minus = 'img_meas_MID00648_FID31367_pulseq_zxminus20.mat';
fn.x2y2.plus = 'img_meas_MID00658_FID31377_pulseq_x2_y2plus20.mat';
fn.x2y2.minus = 'img_meas_MID00656_FID31375_pulseq_x2_y2minus20.mat';
fn.zy.plus = 'img_meas_MID00652_FID31371_pulseq_zyplus20.mat';
fn.zy.minus = 'img_meas_MID00654_FID31373_pulseq_zyminus20.mat';
fn.x.plus = 'img_meas_MID00625_FID31344_pulseq_xplus20.mat';
fn.x.minus = 'img_meas_MID00622_FID31341_pulseq_xminus20.mat';
fn.y.plus = 'img_meas_MID00631_FID31350_pulseq_yplus20.mat';
fn.y.minus = 'img_meas_MID00628_FID31347_pulseq_yminus20.mat';
fn.z.plus = 'img_meas_MID00636_FID31355_pulseq_zplus20.mat';
fn.z.minus = 'img_meas_MID00634_FID31353_pulseq_zminus20.mat';

shims = {'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'};

nx = 60; ny = 60; nz = 50;

% get mask
load(fn.baseline);
img1 = img1((end/2-nx/2):(end/2+nx/2-1), :, :, :, 1);
imsos = sqrt(sum(abs(img1).^2,4));
mask = imsos > 0.1*max(imsos(:));

dte = 1e-3;   % sec
F = zeros(nx,ny,nz,length(shims));
for ii = 1:length(shims)
	load(eval(sprintf('fn.%s.plus', shims{ii})));
	img1 = img1((end/2-nx/2):(end/2+nx/2-1), :, :, :, :);
	pc = phasecontrastmulticoil(img1(:,:,:,:,1), img1(:,:,:,:,3));
	fplus = pc/(2*pi*dte);  % Hz

	load(eval(sprintf('fn.%s.minus', shims{ii})));
	img1 = img1((end/2-nx/2):(end/2+nx/2-1), :, :, :, :);
	pc = phasecontrastmulticoil(img1(:,:,:,:,1), img1(:,:,:,:,3));
	fminus = pc/(2*pi*dte);  % Hz

	f = fplus-fminus;
	f(~mask) = 0;
	F(:,:,:,ii) = f;
	subplot(3,3,ii); im(F(:,:,:,ii)); title(shims{ii}); colormap jet;
	pause(1);
end	

return;


