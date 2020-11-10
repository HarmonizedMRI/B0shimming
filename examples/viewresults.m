load results   % shat, H, A, mask, f0

f = H*A*shat + f0;   % predicted fieldmap after applying shims
f0 = embed(f0,mask);
f = embed(f,mask);
im(cat(1, f0, f)); colorbar;

