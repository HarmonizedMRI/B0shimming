% Create shimvol.mat (define shim region)
% Here we do this by skull stripping using FSL's 'bet' tool

m = magraw;  % magnitude image

if true
    mask = magraw > 0.2*max(magraw(:));
    save shimvol.mat mask
    return
end

% reduce resolution to test bet tool
if false
    d = fftshift(fftn(fftshift(m)));
    n = 92;
    r = (size(d,1)/2-n/2):(size(d,1)/2+n/2-1);
    d = d(r, r, :);
    m = fftshift(ifftn(fftshift(d)));
    m = abs(m);
end

% flip dimensions to match expected orientation for nii images
m = flipdim(m, 1);
m = flipdim(m, 2);
m = flipdim(m, 3);

niftiwrite(m, 'in.nii');

input('do: \`bet in.nii out.nii.gz -f 0.4 -m\`, then press any key to continue');

mb = niftiread('out.nii.gz');
mb = flipdim(mb, 1);
mb = flipdim(mb, 2);
mb = flipdim(mb, 3);
mask = niftiread('out_mask.nii.gz');
mask = flipdim(mask, 1);
mask = flipdim(mask, 2);
mask = flipdim(mask, 3);

save shimvol.mat mask
