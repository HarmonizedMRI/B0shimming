% Create shimvol.mat (define shim region)
% Here we do this by skull stripping using FSL's 'bet' tool

m = abs(x1);  % magnitude image

% flips so the output of bet has the original orientation
m = flipdim(m, 1);
m = flipdim(m, 2);
m = flipdim(m, 3);
niftiwrite(m, 'in.nii');

input('do: \`bet in.nii out.nii.gz -f 0.4 -m\`, then press any key to continue');

mask = niftiread('out_mask.nii.gz');

save shimvol.mat mask
