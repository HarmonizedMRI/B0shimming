function d = loaddata_siemens(data_path)

%% Load data-file obtained with the b0.seq, 

% see writeB0.m.
%
% Input:
%   data_path   [string]    .dat-file name
%
% Output:
%   d    [nx ny nz ncoils 2]    Complex coil data for the 2 echoes

% load data from .dat-file
twix = mapVBVD(data_path);
% twix = data_path;

twix.image.flagRemoveOS  = false; % remove OverSampling
% twix.image.flagRemoveOS  = true; % remove OverSampling


data_unsorted = twix.image.unsorted();
din = permute(data_unsorted, [1 3 2]);   % [nfid nview nslice ncoil]
% twix.image.dataSize()

% discard data during receive gain calibration (see writeB0.m)
% din = din(:,:,2:end,:);  

% """ remove dummy shots """
nzDummy = 1;
nx = size(din,1); 
ny = nx; %twix.hdr.Meas.ImageColumns, twix.hdr.Meas.ImageLines
nz = nx; %#!!!
nCoils = size(din,3);

nTE = 2;   % number of echo times (interleaved)
din = din(:, (nzDummy*ny*nTE+1):end, :);   % [nx ny*nz*nTE nCoils]
%din = reshape(din,[nx,nz,2,ny+1,nCoils]); %???? %
d = zeros(nx, ny*nz, nCoils, nTE);
for iTE = 1:nTE
    d(:, :, :, iTE) = din(:, iTE:nTE:end, :);
end
d = reshape(d, [nx ny nz nCoils nTE]);

% display center kz encode
im(abs(d(:,:,end/2,1,1)).^0.2)
