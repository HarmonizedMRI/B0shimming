
% create shimcal.mat (F S mask FOV)
%makeshimcal;  % edit as needed to ponit to shim calibration Pfiles

% b0init, mask, x1/x2
getb0init;    % edit as needed to point to b0 acquistion

% create f0.mat (f0 FOV)
makef0;       % regularized version of b0init. Uses mask from getb0init.

% create shimvol.mat (mask) by skull stripping x1
makeshimvol;

% Calculate shims
% $ cd ~/github/HarmonizedMRI/B0shimming/julia
% $ julia
% pkg> activate .
% pkg> instantiate
% julia> include("shim.jl");
