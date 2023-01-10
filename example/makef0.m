% Regularize B0 map
%

% regularize field map
% NB! Input to mri_field_map_reg is rad/sec
yik(:,:,:,1) = x1;
yik(:,:,:,2) = x2;
l2b = -1;
fprintf('Regularizing field map...');
load mask  % from getb0init.m
tic
[wmap, wconv] = mri_field_map_reg(yik, [0 dte], ...
    'mask', mask, 'winit', b0init*2*pi, 'l2b', l2b);
fprintf(' done\n');
toc
wmap = wmap .* mask;   % rad/sec

% final b0 estimate
b0 = wmap / (2*pi);  % Hz

f0 = b0;

save f0.mat f0 FOV

