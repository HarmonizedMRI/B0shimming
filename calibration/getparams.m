function p = getparams

%% calibration scan parameters

p.AmpLinear = [-10 10];     % loop over these (integer) hardware amplitudes for each linear shim. Must match shimcal.pl.
p.AmpHO = [-500 500];       % loop over these (integer) hardware amplitudes for each high-order shim.  Must match shimcal.pl.

% NB! If passing 'sys' to writemod.m, 'maxGrad' MUST match the physical system limit -- since gradients are scaled relative to this.
% 'maxSlew' can always be a design choice, i.e., it can be at or below the physical system limit.
% Here we will therefore use 'sys' when designing waveforms, and 'sysGE' when writing them to .mod files with writemod.m.
mxs = 10.0;    % max slew [G/cm/msec]. Go easy to minimize PNS.
p.sys = toppe.systemspecs('maxSlew', mxs, 'slewUnit', 'Gauss/cm/ms', 'maxGrad', 5, 'gradUnit', 'Gauss/cm');

p.res = 0.3;                    % voxel size (cm) (isotropic)
p.fov = 20;                     % in-plane fov (cm) (isotropic)
p.n = 2*round(p.fov/p.res/2);   % matrix size (isotropic)
p.oprbw = 125/4;                % kHz

p.nCyclesSpoil = 1;             % readout spoiler gradient

p.dTE = [0 1 3];                % delta TE (msec)

p.rf.flip = 4;                   % excitation angle (degrees). Low flip for spin density weighting.
p.rf.slabThick = p.fov*0.85;     % cm
p.rf.tbw = 12;                   % time-bandwidth product of SLR pulse 
p.rf.dur = 1;                    % RF pulse duration (msec)
p.rf.ftype = 'min';              % a good option for 3D imaging is 'min'; otherwise 'ls' is common
p.rf.nCyclesSpoil = p.n*p.nCyclesSpoil;      % slice-select spoiler gradient size

%% GUI parameters
p.zeroFillFactor = 4;           % 'zoom' factor for display purposes

