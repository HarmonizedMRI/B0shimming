function reconstructB0(inputData, outputFile, options)
%RECONSTRUCTB0 Prepare B0 reconstruction data for Julia processing.
%
% reconstructB0(scanArchiveFile, outputFile)
% reconstructB0(dCell, outputFile)
% reconstructB0(..., Name, Value)
%
% The first input can be either:
%
%   1. A GE ScanArchive filename:
%
%        reconstructB0('data.h5', 'reconstruction_input.mat')
%
%   2. A cell array of raw acquisitions, such as returned by:
%
%        dCell = pge2.utils.loaddata('data.h5');
%        reconstructB0(dCell, 'reconstruction_input.mat')
%
%      The cell array must contain
%
%        nTE * ny * (NumDummyZ + nz)
%
%      acquisitions. Each cell contains an [nx, nCoils] array.
%
% The output .mat file is used by the Julia reconstruction code to
% calculate an unwrapped and regularized B0 field map.
%
% Name-value options:
%
%   'MatrixSize'       Acquisition matrix size [nx ny nz].
%                      Default: [100 100 100]
%
%   'FOV'              Field of view in meters [x y z].
%                      Default: [0.24 0.24 0.24]
%
%   'NumDummyZ'        Number of dummy z-encodes.
%                      Default: 1
%
%   'MaskThreshold'    Relative magnitude threshold used by Julia.
%                      Default: 0.1
%
%   'B0'               Scanner field strength in Tesla.
%                      Default: 3.0

arguments
    inputData
    outputFile (1,:) char

    options.MatrixSize (1,3) double ...
        {mustBePositive, mustBeInteger} = [100 100 100]

    options.FOV (1,3) double ...
        {mustBePositive} = [0.24 0.24 0.24]

    options.NumDummyZ (1,1) double ...
        {mustBeNonnegative, mustBeInteger} = 1

    options.MaskThreshold (1,1) double ...
        {mustBeInRange(options.MaskThreshold,0,1)} = 0.1

    options.B0 (1,1) double ...
        {mustBePositive} = 3.0
end

% -------------------------------------------------------------------------
% Acquisition parameters
% -------------------------------------------------------------------------

sys = mr.opts('B0', options.B0);

fatChemShift = 3.5e-6;                  % 3.5 ppm
fatOffresFreq = sys.gamma * sys.B0 * fatChemShift;  % Hz

% Fat and water in phase at both echoes.
TE = (1 / fatOffresFreq) * [1 2];       % sec

nTE = length(TE);
assert(nTE == 2, 'Acquisition must have exactly two echoes.');

nx = options.MatrixSize(1);
ny = options.MatrixSize(2);
nz = options.MatrixSize(3);

nDummyZ = options.NumDummyZ;


% -------------------------------------------------------------------------
% Obtain raw acquisitions
% -------------------------------------------------------------------------

if ischar(inputData)

    fprintf('Loading %s...\n', inputData);
    dCell = pge2.utils.loaddata(inputData);

elseif iscell(inputData)

    dCell = inputData;

else

    error(['First input must be either a ScanArchive filename ' ...
           'or a cell array of raw acquisitions.']);

end


% -------------------------------------------------------------------------
% Check raw data
% -------------------------------------------------------------------------

nShots = length(dCell);

expectedShots = nTE * ny * (nDummyZ + nz);

assert( ...
    nShots == expectedShots, ...
    ['The number of TRs in the data (%d) does not match the ' ...
     'expected number (%d).'], ...
    nShots, expectedShots);

assert( ...
    size(dCell{1},1) == nx, ...
    ['MatrixSize(1) = %d does not match the acquired readout ' ...
     'dimension (%d).'], ...
    nx, size(dCell{1},1));

nCoils = size(dCell{1},2);


% -------------------------------------------------------------------------
% Sort acquisitions into k-space array
%
% d has dimensions:
%
%   [nx, nTE*ny, nz, nCoils]
%
% Echoes are interleaved along dimension 2.
% -------------------------------------------------------------------------

d = zeros( ...
    nx, ...
    nCoils, ...
    nTE * ny, ...
    nz, ...
    'like', dCell{1});

% First acquisition after dummy z-encodes.
shot = nDummyZ * nTE * ny + 1;

for z = 1:nz
    for y = 1:nTE*ny
        d(:,:,y,z) = dCell{shot};
        shot = shot + 1;
    end
end

% [nx nCoils nTE*ny nz] -> [nx nTE*ny nz nCoils]
d = permute(d, [1 3 4 2]);


% -------------------------------------------------------------------------
% Save input for Julia processing
% -------------------------------------------------------------------------

fprintf('Saving %s...\n', outputFile);

addpath ../julia/    % saveReconstructionInput.m
saveReconstructionInput( ...
    outputFile, ...
    d, ...
    TE, ...
    'MaskThreshold', options.MaskThreshold);

% Save acquisition geometry for NIfTI output.
fov = single(options.FOV);
matrix_size = int32(options.MatrixSize);

save( ...
    outputFile, ...
    'fov', ...
    'matrix_size', ...
    '-append');

fprintf('Done.\n');

