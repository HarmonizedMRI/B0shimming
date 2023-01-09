% Recon initial b0map, unwrapped

% data file location
[status, tmp] = system('hostname');
hostname = strip(tmp); % remove newline
if strcmp(hostname, 'quickstep')
    datDir = '/mnt/storage/jfnielse/data/20221013_UM3TUHP_3dspiral/';
else
    datDir = '/media/jon/USB/Data/20221013_UM3TUHP_3dspiral/';
end

% fov and matrix size 
res = 0.24;       % cm
N = [92 92 42];   % image matrix size
FOV = N*res;     % cm

% TE shift(s) for b0 mapping
deltaTE = [0 1000/440/2 1000/440];  % TE delays (ms)

pfile = [datDir 'P,b0.7'];
readoutFile = [datDir 'readout_b0.mod'];

% get coil images
echo1 = 1;
echo2 = 3;
[im1, magraw] = toppe.utils.recon3dft(pfile, ...
    'echo', echo1, ...
    'readoutFile', readoutFile, ...
    'alignWithUCS', true);  
im2 = toppe.utils.recon3dft(pfile, ...
    'echo', echo2, ...
    'readoutFile', readoutFile, ...
    'alignWithUCS', true);  

% echo time difference
dte = 1e-3*(deltaTE(echo2) - deltaTE(echo1));  % TE difference, sec

% mask
thr = 0.05;
mask = magraw > thr*max(magraw(:));

mag = magraw.*mask;

% get phase difference map (th)
if false
    if size(im1, 4) > 1   % multicoil
        th = toppe.utils.phasecontrastmulticoil(im2, im1);
    else
        th = angle(im2./im1).*mask;
    end
else
    % complex coil combination
    getsens;
    x1 = coilcombine(im1, sens, mask);
    x2 = coilcombine(im2, sens, mask);
    th = angle(x2./x1).*mask;  % radians. Input to unwrap.
    th(isnan(th)) = 0;
end

b0wrap = th/(2*pi)/dte; % Hz

% unwrap in Julia and save result in unwrap/thuw.mat
cd unwrap
save th th mask mag
cd ..

input('Run unwrap/main.jl, then press any key to continue');
% $ cd unwrap
% $ julia
% Press `]` to enter the Julia package manager.
% pkg> activate .
% pkg> instantiate
% Press backspace to get back to the Julia prompt.
% julia> include("main.jl")    # input to unwrap() is phase image in radians

% load unwrapped phase map
load unwrap/thuw           % radians
thuw = thuw.*mask;         % unwrap() adds pixels at edges!

b0init = thuw/(2*pi)/dte; % Hz

return

% regularize field map in Matlab
% NB! Input to mri_field_map_reg is rad/sec
yik(:,:,:,1) = x1;
yik(:,:,:,2) = x2;
l2b = -1;
fprintf('Regularizing field map...');
tic
[wmap, wconv] = mri_field_map_reg(yik, [0 dte], ...
    'mask', mask, 'winit', thuw/dte, 'l2b', l2b);
fprintf(' done\n');
toc
wmap = wmap .* mask;   % rad/sec

% final b0 estimate
b0 = wmap / (2*pi);  % Hz

return


% old code

if size(im1, 4) > 1   % multicoil
    th = toppe.utils.phasecontrastmulticoil(im2, im1);
else
    th = angle(im2./im1);
end

b0 = th/(2*pi)/(dte*1e-3); % Hz

b0 = b0 .* mask;

save b0.mat b0 mask FOV 
