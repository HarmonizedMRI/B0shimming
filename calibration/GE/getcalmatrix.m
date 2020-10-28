function A = getcalmatrix(fmaps, p, mask)
% function A = getcalmatrix(fmaps, p, mask)
%
% Inputs:
%   fmaps    [nx ny nz 8]    (Hz) Phase-unwrapped 3D fieldmaps obtained by turning on/off each shim term.
%   
%   seq      struct          experimental parameters. See ../getparams.m.
%   mask     [nx ny nz]      ROI to shim over
%
% Output:
%   A        [9 9]   calibration matrix (includes zero-order term)
%
% Test function:
%  >> addpath ..
%  >> seq = getparams();
%  >> A = getCalMatrix('test',seq);

if ischar(fmaps)
	% test this script
	% synthesize noisy fieldmaps
	[nx,ny,nz] = deal(seq.n, seq.n, seq.n);
	[X,Y,Z] = getGrid(nx,ny,nz);
	mask = ones(nx,ny,nz);
	C = getSHbasis(X,Y,Z);
	fmaps = zeros(nx,ny,nz,8);
	for ii = 1:8
		s = zeros(9,1);
		s(ii+1) = 100;
		s(1) = randn(1)/2;
		fmaps(:,:,:,ii) = reshape(C*s, [nx ny nz]);
	end
	fmaps = fmaps + randn(size(fmaps))/40;
end

% image dimensions
[nx,ny,nz,nshim] = size(fmaps);

% grid
[X,Y,Z] = getGrid(nx,ny,nz);

% spherical harmonics basis
C = getSHbasis(X,Y,Z);    % C = [nx*ny*nz 9]

% applied shim values 
% see shimcal.pl. Values are the difference in amplitude between the two acquisitions.
S = zeros(9,9);
S(1,1) = 1;           % zero-order term (center-frequency offset)
for ii = 2:4
	S(ii,ii) = diff(seq.AmpLinear);   % x, y, z
end
for ii = 5:9
	S(ii,ii) = diff(seq.AmpHO);       % z2, xy, zx, x2-y2, zy
end

% calculate calibration matrix
W = spdiag(mask(:));
N = nx*ny*nz;
F = [ones(N,1) reshape(fmaps, [N nshim])];     % [N 9]   Here we add the zero-order (B0) term.
%A = inv(C'*C)*C'*F*inv(S);                    % [9 9]
A = inv(C'*W*W*C)*C'*W*W*F*inv(S);             % [9 9]

return
