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
twix.image.flagRemoveOS  = false; % remove OverSampling
% twix.image.flagRemoveOS  = true; % remove OverSampling


% data = squeeze(twix.image());
data_unsorted = twix.image.unsorted();
din = permute(data_unsorted, [1 3 2]);   % [nfid nview nslice ncoil]
% twix.image.dataSize()


% discard data during receive gain calibration (see writeB0.m)
% din = din(:,:,2:end,:);  
% """ remove dummy shots """
nzDummy = 1; %2; %# should be reduced to 1 in future
nx = size(din,1); 
ny = nx; 
nz = nx; %#attention hard coded shit!!!
nCoils = size(din,3);

% [nfid,ncoil,nslice,necho,nview] = size(din);   % necho = 1 by construction

% nRead = 2*nx; %#attention hard coded shit!!!
% nRead = 2*nx/2; %#attention hard coded shit!!!

% nDummyShots = nzDummy * (nx + ny)+1; %???+1 hm probably not *nCoils
% nDummyShots = nzDummy * ny * nCoils +1; %???+1 hm probably not *nCoils
% data = din(:, nDummyShots:end, :); %# [nFid nCoils 2*ny*nz]
% data = reshape(data,[nx,ny*2,nz,nCoils]); %???? %

% % % data = reshape(data,[nRead,2*ny,nz,2]); %???? %
% % din = data;

% % And here's the k-space for the first coil and the first slice:
% figure,imagesc(abs(squeeze(data(:,:,30,1))).^0.2);

% din = flipdim(din, 1);   % TOPPE data is flipped along first dimension
din = reshape(din,[nx,ny,2*(nz+1),nCoils]); %???? %

% din = permute(din, [1 5 3 2 4]);   % [nfid nview nslice ncoil]
% 
% % discard data during receive gain calibration (see writeB0.m)
din = din(:,:,3:end,:); 
% din = din(:,:,2:end,:); ????

% construct output matrix
% [nx ny nz nCoils] = size(din(:,1:2:end,:,:));
d = zeros(nx, ny, nz, nCoils, 2);
% d(:,:,:,:,1) = din(:,1:2:end,:,:);   % TE1 data, size [60 60 60]???
% d(:,:,:,:,2) = din(:,2:2:end,:,:);   % TE2 data, size [60 60 60]???
d(:,:,:,:,1) = din(:,:,1:2:end,:);   % TE1 data, size [60 60 60]???
d(:,:,:,:,2) = din(:,:,2:2:end,:);   % TE2 data, size [60 60 60]???

% % % And here's the k-space for the first coil, the first slice and TE1:
figure,imagesc(abs(squeeze(d(:,:,30,1,1))).^0.2);
end