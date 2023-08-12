preamble;

if ~exist('b0pre', 'var')
    b0pre = getb0('P,b0.7', 'readout_b0.mod', deltaTE);
end
if ~exist('b0post', 'var')
    b0post = getb0('P,b0,post.7', 'readout_b0.mod', deltaTE);
end

load mask
b0pre = b0pre.*mask;
b0post = b0post.*mask;

[nx ny nz] = size(b0pre);
z = 2:2:N(3);
im(cat(1, b0pre(:,:,z), b0post(:,:,z)), 200*[-1 1]);
colormap jet;
