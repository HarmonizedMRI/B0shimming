function A = getcalmatrix(F, H, S)
% function A = getcalmatrix(F, H, S)
%
% Get 2nd order shim calibration matrix.
% Not that this function 'quietly' inserts the DC (constant B0 offset) terms in F and S
% prior to calculating A.
%
% Inputs:
%   F    [N 8]      (Hz) Phase-unwrapped fieldmaps obtained by turning on/off each of the 8 (3 linear + 5 2nd order) shims.
%                   N = number of voxels ('control points')
%   H    [N 9]      SH basis matrix, evaluated at the same spatial locations (control points) as F.
%                   Unlike F and S, H also contains the DC term (1st column in H). See getSHshims.m.
%   S    [8 8]      (Pairwise difference in) shim currents used to obtain F (hardware units)
%   
% Output:
%   A    [9 9]      calibration matrix (includes DC term), expressed in hardware units

% Add DC term to S and F
Sfull = zeros(9,9);
Sfull(2:end,2:dn) = S;
S = Sfull; 
F = [ones(N,1) F]; 

% get A
A = inv(C'*C)*C'*F*inv(S);                    % [9 9]
%A = inv(C'*W*W*C)*C'*W*W*F*inv(S);             % [9 9]

return
