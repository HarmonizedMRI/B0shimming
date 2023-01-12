% Regularize B0 map

tic
if true
    % regularize field map
    % NB! Input to mri_field_map_reg is rad/sec
    yik(:,:,:,1) = magraw;
    yik(:,:,:,2) = magraw .* exp(1i*2*pi*b0init*dte);
    l2b = -1;
    fprintf('Regularizing field map...');
    load mask  % from getb0init.m
    [wmap, wconv] = mri_field_map_reg(yik, [0 dte], ...
        'mask', mask, 'winit', b0init*2*pi, 'l2b', l2b);
    fprintf(' done\n');
    wmap = wmap .* mask;   % rad/sec

    % final b0 estimate
    b0 = wmap / (2*pi);  % Hz
else
    input('Run b0reg/main.jl, then press any key to continue');
    load b0reg/b0;
end
toc 

f0 = b0;

save f0.mat f0 FOV

