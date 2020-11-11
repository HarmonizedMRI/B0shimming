load results   % shat, H, A, mask, f0

f = H*A*shat + f0;   % predicted fieldmap after applying shims
f0 = embed(f0,mask);
f = embed(f,mask);
im(cat(1, f0, f)); colorbar;

mask2 = mask;
nz = size(mask,3);
mask2(:,:,[1:nz/4 3*nz/4:nz]) = 0;
fs = smooth3(f);
[max(abs(fs(mask))) sqrt(sum(f(mask).^2))]

