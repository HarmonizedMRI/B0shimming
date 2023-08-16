% Create B0 mapping scan files for GE and Siemens

isTest = false;

studyDir= '20221215_UM1.5TSola_shim/';

[status, tmp] = system('hostname');
hostname = strip(tmp); % remove newline
if strcmp(hostname, 'quickstep')
    datDir = '/mnt/storage/jfnielse/data/';   % note trailing slash
else
    datDir = '/media/jon/USB/Data/';
end

addpath ~/github/HarmonizedMRI/Calibration/b0/GE/    % b04ge.m
addpath ~/github/toppeMRI/PulseGEq/                  % +pulsegeq toolbox
%addpath ~/Downloads/Pulseq/pulseq-1.4.0/matlab/      % +mr toolbox
addpath ~/github/pulseq/matlab/      % +mr toolbox

sysGE = toppe.systemspecs('maxSlew', 15, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 152, ... % (us) best if multiple of 4us
    'daqdel', 152, ...  % (us) best if multiple of 4us
    'gradient', 'xrm', ... % xrm: MR750; hrmb: UHP
    'timessi', 200);    % us

% system struct for Siemens
sysSiemens = mr.opts('MaxGrad', sysGE.maxGrad*10, 'GradUnit', 'mT/m', ...
    'MaxSlew', sysGE.maxSlew*10, 'SlewUnit', 'T/m/s', ...
    'B0', 3.0, ...
    'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, ...
    'adcDeadTime', 10e-6);

% Create B0 mapping sequence files for GE (TOPPE)
res = 0.4;  % iso voxel size (cm)
N = [60 60 60];
FOV = N*res;  % cm
flip = 5;             % degrees
deltaTE = [0 1000/440];   % (ms) acquire 2 images with TE extended by 0 and 2.3 ms
