% Create stack of spirals sequence
% ...and calibration (B0) mapping sequence

preamble;

% B0 mapping/coil sensitivity maps
addpath ~/github/HarmonizedMRI/Calibration/b0/GE   % b04ge.m
flip = 4;
b04ge(sys.ge, N, FOV, flip, deltaTE, ...
    'tbw', 8, ...
    'rfDur', 1, ...   % ms
    'ftype', 'min', ...
    'slabThick', FOV(3) * 0.85, ...
    'rfSpoilSeed', 117, ...
    'nCyclesSpoil', 2, ... % spoiling gradient cycles across one voxel
    'fatsat', false);

return

