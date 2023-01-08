% Create the file shimcal.mat needed for B0 shimming.
%
% shimcal.mat   Contains F, S, mask, FOV
% 
% F = [nx ny nz 8] (Hz), in order 'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'
% S = [8 8], shim amplitudes used to obtain F (hardware units)
% mask = [nx ny nz] object support
% FOV = [1 3]  cm
%
% This script reconstructs b0 scans acquired on a GE scanner with TOPPE, 
% by running shimcal.pl. See github/jfnielsen/scanLog/20220907_UM3TMR750_B0shim.
%
% Assumes the readout file is 'readout_b0.mod'

% hardware specs
sys = toppe.systemspecs('maxSlew', 10, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 58, ... % (us) best if multiple of 4us
    'daqdel', 60, ...  % (us) best if multiple of 4us
    'gradient', 'xrm', ... % xrm: MR750; hrmb: UHP
    'timessi', 200);    % us

% FOV and matrix size
FOV = 24 * [1 1 1];  % cm  
N = 60 * [1 1 1]; 

flip = 5;           % excitation angle (degrees)
deltaTE = [0 2.0];  % Change in TE (from minimum) for each scan (ms)

Shims = {'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'};   
AmpLinear = [-10 10];  % see shimcal.pl
AmpHO = [-500 500];    % see shimcal.pl

S = diag([repmat(diff(AmpLinear), [1 3]) repmat(diff(AmpHO), [1 5])]);

clear pfile b0tmp

% get Pfile names
for ii = 1:3
    for jj = 1:length(AmpLinear)
        pfile{ii,jj} = sprintf('P,%s,%d.7', Shims{ii}, AmpLinear(jj));
    end
end

for ii = 4:length(Shims)
    for jj = 1:length(AmpHO)
        pfile{ii,jj} = sprintf('P,%s,%d.7', Shims{ii}, AmpHO(jj));
    end
end

% get matrix size
[~, imsos] = toppe.utils.recon3dft(pfile{1, 1}, ...
    'echo', 1, 'readoutFile', 'readout_b0.mod', 'alignWithUCS', true);
   
F = zeros([size(imsos) size(pfile,1)]);

% get difference fieldmaps for each shim
for ii = 1:size(pfile,1)
    for jj = 1:2
        im1 = toppe.utils.recon3dft(pfile{ii, jj}, ...
            'echo', 1, 'readoutFile', 'readout_b0.mod', 'alignWithUCS', true);
        im2 = toppe.utils.recon3dft(pfile{ii, jj}, ...
            'echo', 2, 'readoutFile', 'readout_b0.mod', 'alignWithUCS', true);
        if size(im1, 4) > 1   % multicoil
            th(:,:,:,jj) = toppe.utils.phasecontrastmulticoil(im2, im1);
        else
            th(:,:,:,jj) = angle(im2./im1);
        end
        b0(:,:,:,jj) = th(:,:,:,jj)/(2*pi)/(diff(deltaTE)*1e-3); % Hz
    end
    F(:,:,:,ii) = diff(b0, 1, 4);
end

% mask
mask = imsos > 0.2*max(imsos(:));

% save to file, to be used by B0shimming Julia code
FOVcal = FOV;
[Xcal, Ycal, Zcal] = toppe.utils.grid2xyz(N, FOV);
save shimcal.mat F S mask FOV
system('rsync -avz shimcal.mat ~/shimtmpfiles/');

