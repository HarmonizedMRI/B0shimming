% Create B0 mapping scan files for GE and Siemens

addpath ~/github/HarmonizedMRI/Calibration/b0/GE/    % b04ge.m
addpath ~/github/toppeMRI/PulseGEq/                  % +pulsegeq toolbox
%addpath ~/Downloads/Pulseq/pulseq-1.4.0/matlab/      % +mr toolbox
addpath ~/github/pulseq/matlab/      % +mr toolbox

params;

if true
	b04ge(sysGE, N, FOV, flip, deltaTE);
	system('cp tipdown.mod tipdown_orig.mod');
else
	system('cp tipdown_orig.mod tipdown.mod');
end

% Convert to Pulseq
% First add non-zero nChop to make Siemens happy
[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod('tipdown.mod');
rf = [rf; zeros(48,1)];
toppe.writemod(sysGE, 'rf', rf, 'gz', gz, 'desc', desc, ...
	'nomflip', flip, ...  % needed to get right B1 scaling in .seq file
    'nChop', [48 48], 'ofname', 'tipdown.mod');
system(sprintf('tar cf b0.tar seqstamp.txt scanloop.txt modules.txt *.mod'));

if false
	pulsegeq.ge2seq('b0.tar', sysGE, sysSiemens, 'seqFile', 'b0.seq');
else
	% test
	pulsegeq.ge2seq('b0.tar', sysGE, sysSiemens, 'seqFile', 'b0.seq', 'nt', 50);
	seq = mr.Sequence(sysSiemens);
	seq.read('b0.seq');
	b = seq.getBlock(1);
	%rep = seq.testReport;
	%fprintf([rep{:}]);
end
