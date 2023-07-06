addpath ~/github/HarmonizedMRI/B0shimming/sequence/Pulseq/
% create .seq file
%writeB0;  

% Convert to TOPPE .tar file for execution on GE scanners
fn = 'b0.seq';

% TOPPE system specs struct
peakRF = max(abs([0.25, rf.signal/sys.gamma*1e4]));  % Gauss
maxGrad = sys.maxGrad/sys.gamma*100;   % G/cm
maxSlew = sys.maxSlew/sys.gamma/10;    % G/cm/ms
sysGE = toppe.systemspecs('maxGrad', maxGrad, ...   % G/cm
    'maxSlew', maxSlew, ...           % G/cm/ms
    'maxRF', peakRF, ...         % Gauss. Must be >= peak RF in sequence.
    'timessi', 100, ...                             % us
    'maxView', 60, ...
    'rfDeadTime', sys.rfDeadTime*1e6, ...           % us
    'rfRingdownTime', sys.rfRingdownTime*1e6, ...   % us
    'adcDeadTime', sys.adcDeadTime*1e6);            % us

verbose = true;    % so it doesn't remove the individual scan files from local folder
pulsegeq.seq2ge(seq, sysGE, 'verbose', verbose, 'tarFile', [fn(1:(end-4)) '.tar']);

return;

% String modules together just like on the scanner, and plot.
% By default, this function reads the scan files in the local folder
% (modules.txt, scanloop.txt, and the .mod files)
nBlocksPerTR = 8;
nTR = 4;
blockFirst = 1; blockLast = nTR*nBlocksPerTR;
figure; 
[rf, gx, gy, gz] = toppe.plotseq(blockFirst, blockLast, sysGE, 'gmax', mxg/10);

% save waveforms
hdr.grad.raster = 4;  % microseconds
hdr.rf.raster = 4;  % microseconds
hdr.grad.unit = 'Gauss/cm';
hdr.rf.unit = 'Gauss';
hdr.script = 'github/jfnielsen/scanLog/20230530_MAGNUS_MariaEngel/diffusion/magnus/main.m';
save([fn(1:(end-4)) '_waveforms.mat'], 'rf', 'gx', 'gy', 'gz', 'hdr');
