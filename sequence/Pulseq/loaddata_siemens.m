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

% twix.image.flagDoAverage = true; %???
twix.image.flagRemoveOS  = true; % remove OverSampling

% data = squeeze(twix.image());
data_unsorted = twix.image.unsorted();
din = permute(data_unsorted, [1 3 2]);   % [nfid nview nslice ncoil]
% twix.image.dataSize()



% discard data during receive gain calibration (see writeB0.m)
% din = din(:,:,2:end,:);  
% """ remove dummy shots """
nzDummy = 2; %# see b04ge.m %#attention hard coded shit!!!    should be reduced to 1 in future
nx = 60;
ny = 60; %#attention hard coded shit!!!
nz = 60;
nCoils = size(din,3);
nRead = 2*nx; %#attention hard coded shit!!!

nDummyShots = nzDummy * ny * nCoils +1; %???+1
data = din(:, nDummyShots:end, :); %# [nFid nCoils 2*ny*nz]
data = reshape(data,[nRead,2*ny,nz,2]); %???? %#attention hard coded shit!!!
din = data;

% % And here's the k-space for the first coil and the first slice:
% figure,imagesc(abs(squeeze(data(:,:,30,1))).^0.2);

% %""" crop fov in x to account for Dwell time being fixed to 4us """
% nc = round(nRead/2);  % center of image
% % print('nc',nc)
% lob = nc-nx/2 +1;
% upb = nc+nx/2;
% din = din(lob:upb, :, :, :);

% construct output matrix
% [nx ny nz nCoils] = size(din(:,1:2:end,:,:));
d = zeros(nx, ny, nz, nCoils, 2);
% d(:,:,:,:,1) = din(:,1:2:end,:,:);   % TE1 data, size [60 60 60]
% d(:,:,:,:,2) = din(:,2:2:end,:,:);   % TE2 data, size [60 60 60]
d(:,:,:,:,1) = din(1:2:end,1:2:end,:,:);   % TE1 data, size [60 60 60]???
d(:,:,:,:,2) = din(2:2:end,2:2:end,:,:);   % TE2 data, size [60 60 60]???

% % And here's the k-space for the first coil, the first slice and TE1:
% figure,imagesc(abs(squeeze(d(:,:,30,1,1))).^0.2);
end