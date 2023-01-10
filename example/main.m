
% create shimcal.mat (F S mask FOV)
%makeshimcal;

% b0init, mask, x1/x2
getb0init;  

% create f0.mat (f0 FOV)
makef0;    

% create shimvol.mat (mask)
makeshimvol;
