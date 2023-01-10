
% create shimcal.mat (F S mask FOV)
makeshimcal;

% create f0.mat (f0 FOV)
getb0init;  % b0init, mask, x1/x2
makef0;    

% create shimvol.mat (mask)
makeshimvol;
