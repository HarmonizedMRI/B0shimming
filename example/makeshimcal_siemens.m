% Create the calibration data file shimcal.mat from Simens raw data.
%
% The data was acquired using the script ???
%
% Assumptions:
%InitialMagnification', 800
%  * 2nd order shim coils
%
%  * balanced acquisitions, i.e., two acquisitions are performed
%    for each shim channel, with amplitudes 'a/2' and '-a/2', respectively.
%
%  * ???P-files are named 'P,<channel>,<amp>.7', where
%      <channel> = 'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', or 'zy'
%      <amp> = shim current amplitude
%
% F = [nx_c ny_c nz_c 8] (Hz), in order 'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'
% S = [8 8], shim amplitudes used to obtain F (hardware units)
% mask_c = [nx_c ny_c nz_c], object support
% FOV_c = [1 3] cm

% location of data files 
% datDir = '~/myDataDir/';
% datDir = '~/myDataDir/shim_test4';
datDir ='~/myDataDir/shim_test_8/'
% datDir = '~/myDataDir/calib_data';

% Acquisition parameters. See ../sequence/Pulseq/writeB0.m.
FOV_c = 24*[1 1 1];  % cm  
nx_c = 60; ny_c = nx_c; nz_c = nx_c;
deltaTE = 2.2369e-3;  % TE difference between the two echoes (sec)
% deltaTE = 1000/440 *1e-3; % NW from python script.. !! Attention hard coded shit! 

% Shim channel names and amplitude settings
shims = {'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'};   
AmpLinear = [-10 10];  % see shimcal??
AmpHO = [-100 100];    % see shimcal??
nShim = length(shims);

S = diag([repmat(diff(AmpLinear), [1 3]) repmat(diff(AmpHO), [1 5])]);

%% get difference fieldmaps for each shim (and mask)
F = zeros([nx_c ny_c nz_c nShim]);
nF = 0;
pattern='*.dat';
D=dir([datDir filesep pattern]);
[~,I]=sort(string({D(:).name}));

% b0 = zeros();
for ii = 1:nShim
    for jj = 1:2
        nF = nF + 1;
        data_file_path=[datDir filesep D(I(nF)).name];
        d = loaddata_siemens(data_file_path);
        [b0(:,:,:,jj), mask_c] = reconB0(d, deltaTE, 0.1);
    end
    F(:,:,:,ii) = diff(b0, 1, 4);
end

% % save to file
save shimcal.mat F S FOV_c mask_c
