% Recon initial b0map, unwrapped
% Also save the following that are used by b0reg/main.jl:
% images.mat     complex coil images
% echotime.mat   echo times (sec)
% 
%

if false
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
else
    preamble;
end

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

echotime = deltaTE([echo1 echo2])*1e-3;  % sec
save echotime echotime

clear images;
images(:,:,:,:,1) = im1;
images(:,:,:,:,2) = im2;
save images images

% echo time difference
dte = 1e-3*(deltaTE(echo2) - deltaTE(echo1));  % TE difference, sec

% mask
thr = 0.05;
mask = magraw > thr*max(magraw(:));
save mask mask

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

save b0init b0init

