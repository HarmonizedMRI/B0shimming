function reconstructB0(scanArchiveFile, outputFile, options)
%RECONSTRUCTB0 Load B0 ScanArchive data and prepare input for Julia.
%
% reconstructB0(scanArchiveFile, outputFile)
% reconstructB0(..., Name=Value)
%
% Loads a GE ScanArchive data file acquired with b0.seq, sorts the raw
% k-space data, and saves a MATLAB file for subsequent B0 field-map
% estimation in ../julia/.
%
% Example:
%
%   reconstructB0("data.h5", "reconstruction_input.mat")
%
% Name-value options:
%   MatrixSize      Acquisition matrix size [nx ny nz].
%                   Default: [100 100 100]
%
%   FOV             Field of view in meters [x y z].
%                   Default: [0.24 0.24 0.24]
%
%   NumDummyZ       Number of dummy z-encodes.
%                   Default: 1
%
%   MaskThreshold   Relative magnitude threshold used by the Julia
%                   reconstruction.
%                   Default: 0.1
%
%   B0              Scanner field strength in Tesla.
%                   Default: 3.0

arguments
    scanArchiveFile (1,:) char
    outputFile      (1,:) char

    options.MatrixSize (1,3) double {mustBePositive, mustBeInteger} = [100 100 100]
    options.FOV (1,3) double {mustBePositive} = [0.24 0.24 0.24]
    options.NumDummyZ (1,1) double {mustBeNonnegative, mustBeInteger} = 1
    options.MaskThreshold (1,1) double {mustBeInRange(options.MaskThreshold,0,1)} = 0.1
    options.B0 (1,1) double {mustBePositive} = 3.0
end

% -------------------------------------------------------------------------
% Acquisition parameters
% -------------------------------------------------------------------------

sys = mr.opts('B0', options.B0);

fatChemShift = 3.5e-6;   % 3.5 ppm
fatOffresFreq = sys.gamma * sys.B0 * fatChemShift;  % Hz

% Choose echo times so fat and water are in phase at both echoes.
TE = (1 / fatOffresFreq) * [1 2];   % seconds

nTE = length(TE);
assert(nTE == 2, 'Acquisition must have exactly two echoes.');

nx = options.MatrixSize(1);
ny = options.MatrixSize(2);
nz = options.MatrixSize(3);

nDummyZ = options.NumDummyZ;

% -------------------------------------------------------------------------
% Load GE ScanArchive data
% -------------------------------------------------------------------------

fprintf('Loading %s...\n', scanArchiveFile);

dCell = pge2.utils.loaddata(scanArchiveFile);

nShots = length(dCell);

expectedShots = nTE * ny * (nDummyZ + nz);

assert( ...
    nShots == expectedShots, ...
    ['The number of TRs in the data file (%d) does not match the ' ...
     'expected number (%d).'], ...
    nShots, expectedShots);

assert( ...
    nx == size(dCell{1}, 1), ...
    ['MatrixSize(1) = %d does not match the acquired readout ' ...
     'dimension (%d).'], ...
    nx, size(dCell{1}, 1));

nCoils = size(dCell{1}, 2);

% -------------------------------------------------------------------------
% Sort acquisitions into k-space array
%
% Output organization:
%
%   d(nx, nTE*ny, nz, nCoils)
%
% Echoes are interleaved along dimension 2.
% -------------------------------------------------------------------------

d = zeros(nx, nCoils, nTE * ny, nz, 'like', dCell{1});

shot = nDummyZ * nTE * ny + 1;

for z = 1:nz
    for y = 1:nTE * ny
        d(:, :, y, z) = dCell{shot};
        shot = shot + 1;
    end
end

d = permute(d, [1 3 4 2]);

% -------------------------------------------------------------------------
% Save input for Julia processing
% -------------------------------------------------------------------------

fprintf('Saving %s...\n', outputFile);

saveReconstructionInput( ...
    outputFile, ...
    d, ...
    TE, ...
    MaskThreshold=options.MaskThreshold);

fprintf('Done.\n');


