d = loaddata_siemens(data_path)
% d = loaddata_ge(pfile)
%
% Load GE P-file obtained with the TOPPE version of b0.seq, 
% see writeB0.m. See also the 'Pulseq on GE' manual.
%
% Input:
%   pfile   [string]    GE P-file name
%
% Output:
%   d    [nx ny nz ncoils 2]    Complex coil data for the 2 echoes

% % load data from P-file (set CV2 = 1 to save P-file)
% din = toppe.utils.loadpfile(pfile);
% din = flipdim(din, 1);   % TOPPE data is flipped along first dimension
% [nfid,ncoil,nslice,necho,nview] = size(d);   % necho = 1 by construction
% din = permute(din, [1 5 3 2 4]);   % [nfid nview nslice ncoil]

%% load data from .dat-file
% clear all
% twix = mapVBVD('/home/wehkamp/git/shim_test/2023-02-15-201930.dat');
twix = mapVBVD(data_path);

% twix.image.flagDoAverage = true; %???
twix.image.flagRemoveOS  = true; %???
% twix.image.flagIgnoreSeg = true;
% twix.image.squeeze = true;
% data = squeeze(twix.image());
data_unsorted = twix.image.unsorted();
din = permute(data_unsorted, [1 3 2]);   % [nfid nview nslice ncoil]
% 
% twix.image.dataSize()

% din = flip(data, 1); %NW does this do anything???  % TOPPE data is flipped along first dimension
% [nfid,ncoil,nslice,necho,nview] = size(d);   % necho = 1 by construction

% %                 Col Cha Lin Par Sli Ave Phs Eco Rep Set Seg 
% data = twix.image(  :,  :,  :,  1,  1,  1,  1,  1,  1,  1, :);

%% discard data during receive gain calibration (see writeB0.m)
% din = din(:,:,2:end,:);  
% """ remove dummy shots """
nEcho = 2; %#length(deltaTE) #attention hard coded shit!!!:w
nzDummy = 2; %# see b04ge.m
nx = 60;
ny = 60; %#attention hard coded shit!!!
nz = 60;
nCoils = 2; %#attention hard coded shit!!!
nRead = 2*nx;

nDummyShots = nzDummy * ny * nCoils +1; %???+1
data = din(:, nDummyShots:end, :); %# [nFid nCoils 2*ny*nz]
data = reshape(data,[120,120,60,2]); %???? %#attention hard coded shit!!!
din = data;

%""" crop fov in x to account for Dwell time being fixed to 4us """
nc = round(nRead/2);  %# center of image
% print('nc',nc)
lob = nc-nx/2 +1;
upb = nc+nx/2;
din = din(lob:upb, :, :, :);




%% construct output matrix
% [nx ny nz nCoils] = size(din(:,1:2:end,:,:));
d = zeros(nx, ny, nz, nCoils, 2);
d(:,:,:,:,1) = din(:,1:2:end,:,:);   % TE1 data, size [60 60 60]
d(:,:,:,:,2) = din(:,2:2:end,:,:);   % TE2 data, size [60 60 60]

%% And here's the k-space for the first coil and the first slice:
figure,imagesc(abs(squeeze(d(:,1,:,1,2))).^0.2);