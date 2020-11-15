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
%   H    [N ...]    SH basis matrix, evaluated at the same spatial locations (control points) as F.
%                   Unlike F and S, H also contains the DC term (1st column in H). See getSHshims.m.
%   S    [8 8]      (Pairwise difference in) shim currents used to obtain F (hardware units)
%   
% Output:
%   A    [9 9]      calibration matrix (includes DC term), expressed in hardware units

N = size(F,1);

% Add DC offset term (column) to F and S
% We do it here so the user doesn't have to remember to do it.
F = [ones(N,1) F]; 
Sin = S;
S = zeros(9);
S(1,1) = 1;
S(2:end,2:end) = Sin;

% get A
A = inv(H'*H)*H'*F*inv(S);                    % [9 9]
%A = inv(C'*W*W*C)*C'*W*W*F*inv(S);             % [9 9]

return
