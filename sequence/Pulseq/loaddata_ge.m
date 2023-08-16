function d = loaddata_ge(pfile)
% function d = loaddata_ge(pfile)
%
% Load GE P-file obtained with the TOPPE version of b0.seq, 
% see writeB0.m. See also the 'Pulseq on GE' manual.
%
% Input:
%   pfile   [string]    GE P-file name
%
% Output:
%   d    [nx ny nz ncoils 2]    Complex coil data for the 2 echoes

% load data from P-file (set CV2 = 1 to save P-file)
din = toppe.utils.loadpfile(pfile);
din = flipdim(din, 1);   % TOPPE data is flipped along first dimension
[nfid,ncoil,nslice,necho,nview] = size(din);   % necho = 1 by construction
din = permute(din, [1 5 3 2 4]);   % [nfid nview nslice ncoil]

% discard data during receive gain calibration (see writeB0.m)
din = din(:,:,2:end,:);  

% construct output matrix
[nx ny nz ncoils] = size(din(:,1:2:end,:,:));
d = zeros(nx, ny, nz, ncoils, 2);
d(:,:,:,:,1) = din(:,1:2:end,:,:);   % TE1 data, size [60 60 60]
d(:,:,:,:,2) = din(:,2:2:end,:,:);   % TE2 data, size [60 60 60]
