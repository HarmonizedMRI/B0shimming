function saveReconstructionInput(filename, d, TE, options)
%SAVERECONSTRUCTIONINPUT Save B0 reconstruction input for Julia.
%
% saveReconstructionInput(filename, d, TE)
% saveReconstructionInput(..., Mask=mask, MaskThreshold=0.4)
%
% Inputs
%   filename   Output .mat filename.
%   d          Complex k-space [nx, nEcho*ny, nz, ncoil].
%   TE         Echo times in seconds.
%
% Name-value options
%   Mask           Optional 3D processing mask.
%   MaskThreshold  Relative magnitude threshold used by Julia when Mask is
%                  omitted. Default: 0.4.

arguments
    filename (1,1) string
    d {mustBeNumeric}
    TE double
    options.Mask = []
    options.MaskThreshold (1,1) double {mustBeInRange(options.MaskThreshold,0,1)} = 0.1
end

assert(ndims(d) == 4, ...
    'd must have dimensions [nx, nEcho*ny, nz, ncoil].');
assert(length(TE) == 2, ...
    'The current Julia reconstruction requires exactly two echoes.');
assert(mod(size(d,2), length(TE)) == 0, ...
    'The second dimension of d must be divisible by the number of echoes.');

kspace = single(d);
echo_times = single(TE(:));
mask_threshold = single(options.MaskThreshold);

if isempty(options.Mask)
    save(filename, 'kspace', 'echo_times', 'mask_threshold', '-v7.3');
else
    mask = logical(options.Mask);
    expectedSize = [size(d,1), size(d,2)/length(TE), size(d,3)];
    assert(isequal(size(mask), expectedSize), ...
        'Mask size must be [%d %d %d].', expectedSize);
    save(filename, 'kspace', 'echo_times', 'mask_threshold', 'mask', '-v7.3');
end
end
