% Create matched Pulseq and TOPPE scans for 3D B0 field mapping
%
% Outputs:
%  Pulseq scan file        B0scan.seq
%  TOPPE scan files        modules.txt, scanloop.txt, tipdown.mod, and readout.mod
%
% Usage:
% >> makescanfiles;

% Set paths to Pulseq and TOPPE libraries
%addpath ~/github/pulseq/matlab/         % +mr package
%addpath ~/github/toppeMRI/toppe/        % +toppe package
%addpath ~/github/toppeMRI/PulseGEq/     % +pulsegeq package (Pulseq <--> TOPPE conversion)
% paths for MZ
%addpath ~/pulseq_home/github/pulseq/matlab/
%addpath ~/pulseq_home/github/toppe/
%addpath ~/pulseq_home/github/PulseGEq/

% Settings for GE
GEfilePath = '/usr/g/research/pulseq/';    % Scan file path (GE only)
gmaxGE = 8;                                % Physical hardware spec (Gauss/cm)


%% Acquisition parameters
% Minimum TR will be calculated below
nx = 60;
ny = 60;
fov = [24 24 20];                  % cm
if fov(1) ~= fov(2)
	error('In-plane fov must be square');
end
nz = 2*round(nx*fov(3)/fov(1)/2);  % isotropic voxels
deltaTE = [0 1.0];
deltaTE = [0 0.5 1.0 2.0];       % Change in TE (from minimum) for each of the >= 2 scans needed to estimate B0 field (msec)
oprbw = 125/4;           % Acquisition bandwidth (kHz) for TOPPE scan. Determines gradient readout trapezoid shape.
nCycleSpoil = 2;         % readout spoiler gradient area (cycles across voxel dimension)
ex.flip = 5;             % degrees
ex.thick = 0.85*fov(3);  % excited slab thickness. To avoid wrap-around in z.
ex.tbw = 8;              % time-bandwidth product of excitation pulse
ex.dur = 2;              % duration (msec)

% Set hardware limits used for waveform DESIGN 
limits.design = toppe.systemspecs('maxSlew', 10, 'slewUnit', 'Gauss/cm/ms', ...
	'maxGrad', 2.8, 'gradUnit', 'Gauss/cm', ...
	'maxRf', 0.25, 'rfUnit', 'Gauss');

% Define the PHYSICAL limits for each scanner
% NB! When creating .mod files with toppe.writemod, 'maxGrad' MUST match the PHYSICAL system limit since gradients are scaled relative to this.
ge.system = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
	'maxGrad', gmaxGE, 'gradUnit', 'Gauss/cm', ...
	'maxRf', 0.25, 'rfUnit', 'Gauss');

siemens.system = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', ...
	 'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

%% Create waveforms for TOPPE and associated files (modules.txt, *.mod files). In TOPPE, only 'arbitrary' waveforms are used.

% Create temporary modules.txt (must do it here since this file is used for intermediate steps below)
modFileText = ['' ...
'Total number of unique cores\n' ...
'2\n' ...
'fname	duration(us)	hasRF?	hasDAQ?\n' ...
'readout.mod	0	0	1\n' ...                        % Entries are tab-separated
'tipdown.mod	0	1	0' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);

% Design rf waveform and write to .mod file
ex.nCycleSpoil = nCycleSpoil * ex.thick/(fov(3)/nz);
[ex.rf, ex.g] = toppe.utils.rf.makeslr(ex.flip, ex.thick, ex.tbw, ex.dur, ex.nCycleSpoil, ...
	'type', 'st', ...        % 'st' = small-tip. 'ex' = 90 degree design
	'ftype', 'ls', ...      % minimum-phase SLR design is well suited for 3D slab excitation (can tolerate mild quadratic phase along z)
	'system', limits.design, ...
	'writeModFile', false);
ex.rf = toppe.utils.makeGElength(ex.rf);
ex.g = toppe.utils.makeGElength(ex.g);
toppe.writemod('rf', ex.rf, 'gz', ex.g, ...
	'ofname', 'tipdown.mod', ...
	'system', ge.system);

% Design readout gradients and write to .mod file
% Remember: makegre() creates y phase-encodes based on isotropic in-plane resolution
dz = fov(3)/nz;   % z voxel size (cm)
[acq.gx, acq.gy, acq.gz] = toppe.utils.makegre(fov(1), nx, dz, ...
	'system', limits.design, 'ncycles', nCycleSpoil, ...
	'slewDerate', 0.7, ...   % reduce PNS
	'ofname', 'tmp.mod', ...
	'oprbw', oprbw);
[rf,gx,gy,gz,desc,paramsint16,paramsfloat] = toppe.readmod('tmp.mod');
toppe.writemod('gx', acq.gx, 'gy', acq.gy, 'gz', acq.gz, 'ofname', 'readout.mod', ...
	'system', ge.system, ...
	'desc', desc, 'hdrints', paramsint16, 'hdrfloats', oprbw);

% Determine TR
% We do this by writing a short scanloop.txt file and calculating the resulting minimum TR from that.
toppe.write2loop('setup', 'version', 3);
toppe.write2loop('tipdown.mod', 'textra', max(deltaTE));
toppe.write2loop('readout.mod', 'textra', min(deltaTE));
toppe.write2loop('finish');
% At this point you can view one TR of the TOPPE sequence by typing 'toppe.plotseq(1,2);'
TR = toppe.getTRtime(1,2)*1e3;      % msec


%% Prepare Pulseq sequence object and waveforms

% Create a new Pulseq sequence object
seq = mr.Sequence(siemens.system);         

% Create waveforms suitable for Pulseq by converting units and interpolating.
% Discard zeros at beginning and end of rf waveform.
I = find(abs(ex.rf) == 0);
iStart = find(diff(I)>1) + 1;
iStop = I(iStart+1);
siemens.ex.rf = pulsegeq.rf2pulseq(ex.rf(iStart:iStop), ge.system.raster, seq);  % Gauss -> Hz; 4us -> 1us.
siemens.ex.rfdelay = roundToRaster(iStart * ge.system.raster, siemens.system.gradRasterTime);
siemens.ex.gdelay = max(0, siemens.system.rfDeadTime - siemens.ex.rfdelay);
siemens.ex.g  = pulsegeq.g2pulseq(ex.g, ge.system.raster, seq);    % Gauss/cm -> Hz/m; 4us -> 10us
siemens.acq.gx = pulsegeq.g2pulseq(acq.gx, ge.system.raster, seq); % readout gradient
siemens.acq.gy = pulsegeq.g2pulseq(acq.gy, ge.system.raster, seq); % phase-encode gradient
siemens.acq.gz = pulsegeq.g2pulseq(acq.gz, ge.system.raster, seq); % partition-encode gradient

tmp.rf = mr.makeArbitraryRf(siemens.ex.rf, ex.flip/180*pi, 'system', siemens.system, 'delay', siemens.ex.rfdelay);
tmp.readout = mr.makeArbitraryGrad('z', siemens.acq.gx, siemens.system);
tr_min = mr.calcDuration(tmp.rf) + mr.calcDuration(tmp.readout);   % sec
siemens.delay = roundToRaster(TR*1e-3-tr_min, siemens.system.gradRasterTime);       % sec

% Create Pulseq adc object
[~,~,~,~,~,paramsint16] = toppe.readmod('readout.mod');
nPre = paramsint16(1);  % number of samples in gradient pre-winder and ramp to plateau
nPlateau = paramsint16(2);
acq.preDelay = nPre*ge.system.raster;            % sec
acq.flat = nPlateau*ge.system.raster;            % duration of flat portion of readout (sec)
siemens.acq.N = 2*round(acq.flat/siemens.system.gradRasterTime/2);   % number of readout samples (Siemens)
siemens.acq.dur = siemens.acq.N*siemens.system.gradRasterTime;
pulseq.adc = mr.makeAdc(siemens.acq.N, 'Duration', siemens.acq.dur, 'Delay', acq.preDelay, 'system', siemens.system);

% Create other Pulseq objects that don't need updating in scan loop (except phase)
pulseq.acq.gx = mr.makeArbitraryGrad('x', siemens.acq.gx, siemens.system); 
pulseq.ex.rf = mr.makeArbitraryRf(siemens.ex.rf, ex.flip/180*pi, ...
	'system', siemens.system, 'delay', siemens.ex.rfdelay);
pulseq.ex.g  = mr.makeArbitraryGrad('z', siemens.ex.g, siemens.system, 'delay', siemens.ex.gdelay);


%% Scan loop. In each iteration, 'scanloop.txt' is updated, and blocks are added to seq object
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

toppe.write2loop('setup', 'version', 3);   % Initialize scanloop.txt
for iz = 0:nz     % iz = 0 is used as discarded acquisitions to reach steady state    
	for ii = 1:30
		fprintf('\b');
	end
	fprintf('Writing z-encode %d of %d', iz, nz);
	if iz > 0   
		dabmode = 'on';
	else
		dabmode = 'off';
	end

	% phase-encode loop
	for iy = 1:ny 

		% echo-time loop
		for ie = 1:length(deltaTE)

			% y/z phase-encode amplitudes, scaled to (-1,1)
			a_gy = (iz > 0) * ((iy-1+0.5)-ny/2)/(ny/2); % * ny/nx; 
			a_gz = (iz > 0) * ((iz-1+0.5)-nz/2)/(nz/2);

			% rf excitation (TOPPE)
	  		toppe.write2loop('tipdown.mod', 'RFamplitude', 1.0, 'RFphase', rfphs, ...
				'textra', deltaTE(ie));            % dead time at end of module (msec)

			% rf excitation (Pulseq)
			pulseq.ex.rf.phaseOffset = rfphs;
			seq.addBlock(pulseq.ex.rf, pulseq.ex.g);
			if deltaTE(ie) > 0
				seq.addBlock(mr.makeDelay(1e-3*deltaTE(ie)));
			end

		 	% readout (TOPPE)
			toppe.write2loop('readout.mod', ...
				'Gamplitude', [1 a_gy a_gz]', ...
				'DAQphase', rfphs, ...
				'slice', max(iz,1), 'echo', ie, 'view', iy, ...  % data is stored in 'slice', 'echo', and 'view' indeces. Will change to ScanArchive in future.
				'textra', max(deltaTE)-deltaTE(ie), ...
				'dabmode', dabmode);

			% readout (Pulseq)
			gy = mr.makeArbitraryGrad('y', a_gy*siemens.acq.gy, siemens.system); 
			gz = mr.makeArbitraryGrad('z', a_gz*siemens.acq.gz, siemens.system); 
			pulseq.adc.phaseOffset = rfphs;   % radians
			seq.addBlock(pulseq.acq.gx, gy, gz, pulseq.adc); 
			if deltaTE(ie) < max(deltaTE)
				seq.addBlock(mr.makeDelay(1e-3*(max(deltaTE)-deltaTE(ie))));
			end

			% update rf phase (RF spoiling)
			rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
			rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
		end
	end
end
toppe.write2loop('finish');
fprintf('\n');


%% Write Pulseq file

% check whether the timing of the sequence is correct
fprintf('Checking Pulseq timing...');
[ok, error_report]=seq.checkTiming;
fprintf('\n');

if (ok)
    fprintf('\tTiming check passed successfully\n');
else
    fprintf('\tTiming check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% Write to Pulseq file
fprintf('Writing Pulseq file...');
seq.setDefinition('FOV', fov*1e-2);   % m
seq.setDefinition('Name', '3D B0 mapping');
seq.write('B0scan.seq');
fprintf('done\n');
% parsemr('B0scan.seq');

if false
% Display (part of) sequences
fprintf('Displaying sequences...');
nModsPerTR = 2;    % number of TOPPE modules per TR
nTR = ny/2;        % number of TRs to display
nStart = nModsPerTR * floor(ny*length(deltaTE)+ny/2-nTR/2);
toppe.plotseq(nStart, nStart + nTR*nModsPerTR);
tStart = nStart/nModsPerTR*TR*1e-3;    % sec
tStop = tStart + nTR*TR*1e-3;          % sec
seq.plot('timeRange', [tStart tStop]);
fprintf('\n');

% Display TOPPE sequence in loop/movie mode
%figure; toppe.playseq(nModsPerTR, 'drawpause', false);  
end


%% Create modules.txt and toppe0.meta and write TOPPE files to a tar file

% Create (new) modules.txt 
modFileText = ['' ...
'Total number of unique cores\n' ...
'2\n' ...
'fname	duration(us)	hasRF?	hasDAQ?\n' ...
GEfilePath 'readout.mod	0	0	1\n' ...                        % Entries are tab-separated
GEfilePath 'tipdown.mod	0	1	0' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);

% Create toppe0.meta. This file is the main entry point for toppev3, 
% and must be placed in /usr/g/bin/ (or another hardcoded path). 
metaFileText = ['' ...
GEfilePath 'modules.txt\n' ...
GEfilePath 'scanloop.txt\n' ...
GEfilePath 'tipdown.mod\n' ...
GEfilePath 'readout.mod\n' ];
fid = fopen('toppe0.meta', 'wt');
fprintf(fid, metaFileText);
fclose(fid);

system(sprintf('tar czf B0scan_gmax%d.tgz modules.txt scanloop.txt tipdown.mod readout.mod toppe0.meta', gmaxGE*10));

instr = ['\nPlace toppe0.meta on scanner host (path is hardcoded in the binary)\n' ...
'Untar B0scan.tgz in ' GEfilePath ' on scanner host\n' ];
fprintf(instr);
