% Create the calibration data file shimcal.mat from Simens raw data.
%
% The data was acquired using the script ???
%
% Assumptions:
%
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

% Development notes: See also github/jfnielsen/scanLog/20220907_UM3TMR750_B0shim.

% location of data files 
% datDir = '~/myDataDir/';
datDir = '/home/wehkamp/myDataDir/shim_test';

% Acquisition parameters. See ../sequence/Pulseq/writeB0.m.
FOV_c = 24*[1 1 1];  % cm  
nx_c = 60; ny_c = nx_c; nz_c = nx_c;
% deltaTE = 2.2369e-3;  % TE difference between the two echoes (sec)
deltaTE = 1000/440 *1e-3; % NW from python script.. !! Attention hard coded shit! 

% Shim channel names and amplitude settings
shims = {'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'};   
% AmpLinear = [-10 10];  % see shimcal.pl   %NW GE
% AmpHO = [-500 500];    % see shimcal.pl   %NW GE
AmpLinear = [-20 20];  % see shimcal??   !!!Attention change to 10 ??
AmpHO = [-200 200];    % see shimcal??     !!!Attention change to 100 ???
nShim = length(shims);

S = diag([repmat(diff(AmpLinear), [1 3]) repmat(diff(AmpHO), [1 5])]);

%% test
% file_name = '/2023-02-15-201930.dat';
% data_file_path = append(datDir, file_name);
% .dat-file names (Assuming they were measured in order)

%% Load data from the given directory sorted by name
% path='../data'; % directory to be scanned for data files
% nF=1; % the number of the file / data set to load

% pattern='*.dat';
% D=dir([datDir filesep pattern]);
% [~,I]=sort(string({D(:).name}));
% data_file_path=[datDir filesep D(I(nF)).name];

% fprintf(['loading `' data_file_path '´ ...\n']);
% twix_obj = mapVBVD(data_file_path);

%% get difference fieldmaps for each shim (and mask)
F = zeros([nx_c ny_c nz_c nShim]);
nF = 0;
pattern='*.dat';
D=dir([datDir filesep pattern]);
[~,I]=sort(string({D(:).name}));

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


%% Load the latest fildat_namee from a dir
%path='../IceNIH_RawSend/'; % directory to be scanned for data files
%pattern='*.dat';
%D=dir([path filesep pattern]);
%[~,I]=sort([D(:).datenum]);
%data_file_path=[path D(I(end-5)).name]; % use end-1 to reconstruct the second-last data set, etc.

%% P-file names
% for ii = 1:3
%     for jj = 1:length(AmpLinear)
%         pfile{ii,jj} = sprintf('%s/P,%s,%d.7', datDir, shims{ii}, AmpLinear(jj));
%     end
% end
% 
% for ii = 4:nShim
%     for jj = 1:length(AmpHO)
%         pfile{ii,jj} = sprintf('%s/P,%s,%d.7', datDir, shims{ii}, AmpHO(jj));
%     end
% end