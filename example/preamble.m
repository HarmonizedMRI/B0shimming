% Common protocol settings

[status, tmp] = system('hostname');
hostname = strip(tmp); % remove newline
if strcmp(hostname, 'quickstep')
    datDir = '/mnt/storage/jfnielse/data/20221013_UM3TUHP_3dspiral/';
else
    datDir = '/media/jon/USB/Data/20221013_UM3TUHP_3dspiral/';
    datDir = '';
end

% system hardware limits
sys.ge = toppe.systemspecs('maxSlew', 15, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 152, ... % (us) best if multiple of 4us
    'daqdel', 152, ...  % (us) best if multiple of 4us
    'gradient', 'hrmb', ... % xrm: MR750; hrmb: UHP; hrmw: Premier
    'timessi', 200);    % us

%sys.siemens = mr.opts('MaxSlew', 150, 'SlewUnit', 'T/m/s', ... 
%    'MaxGrad', 28, 'GradUnit', 'mT/m', ...
%    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% fov and matrix size 
FOV = [24 24 24];
N = [60 60 60];   % image matrix size

% TE shift(s) for b0 mapping
deltaTE = [0 1000/440];  % TE delays (ms)


