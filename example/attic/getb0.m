function b0 = getb0(pfile, readoutFile, deltaTE)
%
% Reconstruct b0 map from pfile, first two echoes.
%
%  pfile    string   Pfile name
%  readout  string   e.g., 'readout.mod'
%  deltaTE  [1 2]    echo times (ms)

% get coil images
echo1 = 1;
echo2 = 2;
[im1, magraw] = toppe.utils.recon3dft(pfile, ...
    'echo', echo1, ...
    'readoutFile', readoutFile, ...
    'alignWithUCS', true);  
im2 = toppe.utils.recon3dft(pfile, ...
    'echo', echo2, ...
    'readoutFile', readoutFile, ...
    'alignWithUCS', true);  

dte = 1e-3*(deltaTE(echo2) - deltaTE(echo1));  % TE difference, sec

% mask
thr = 0.05;
mask = magraw > thr*max(magraw(:));
save mask mask

mag = magraw.*mask;

% get phase difference map (th)
if size(im1, 4) > 1   % multicoil
    th = toppe.utils.phasecontrastmulticoil(im2, im1);
else
    th = angle(im2./im1).*mask;
end

b0 = th/(2*pi)/dte; % Hz

