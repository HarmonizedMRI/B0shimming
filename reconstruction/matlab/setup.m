% mr.opts()
system('git clone --branch master git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

% pge2.utils.loaddata() for loading GE ScanArchive files
system('git clone git@github.com:HarmonizedMRI/pge2.git');
addpath pge2/matlab

% ift3.m
system('git clone git@github.com:HarmonizedMRI/utils.git');
addpath utils

% To load the ScanArchive raw data files you will also need the Orchestra toolbox
% which is available for download at http://weconnect.gehealthcare.com/
addpath ~/Programs/orchestra-sdk-2.1-1.matlab/

% We will need textprogressbar.m
system('git clone git@github.com:HarmonizedMRI/pulseg.git');
addpath pulseg/matlab
addpath(genpath('pulseg/matlab/third_party'));

% get 'im' display function (optional)
system('git clone --depth 1 git@github.com:JeffFessler/mirt.git');
cd mirt; setup; cd ..;
