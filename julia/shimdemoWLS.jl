using JLD2
using MIRT: embed!, jim
using Random
using LinearAlgebra

include("getSHbasis.jl")
include("getcalmatrix.jl")

# load calibration data
@load "CalibrationData.jld2" F S fov meta 

(nx,ny,nz,nShim) = size(F)

(x,y,z) = ndgrid(
	range(-fov[1], fov[1], length=nx), 
	range(-fov[2], fov[2], length=ny), 
	range(-fov[3], fov[3], length=nz) 
	)

# create mask
fm = sum(abs.(F), dims=4)[:,:,:,1]
mask = fm .> 200    # as always in Julia, note the '.' (broadcasting)
N = sum(mask[:])

# mask and reshape to [N 8] 
Fm = zeros(N, nShim)    # m for 'masked'
for ii = 1:nShim
	f1 = F[:,:,:,ii]
	Fm[:,ii] = f1[mask]
end

# Get spherical harmonic basis of degree l
l = 2     
H = getSHbasis(x[mask], y[mask], z[mask], l)   # size is [N sum(2*(0:l) .+ 1)]

# Get calibration matrix A
A = getcalmatrix(Fm, H, S);

# Synthesize an example B0 map and optimize shims (minimized RMS residual) for that fieldmap.
shimNum = 5;
f0 = F[:,:,:,4] + 0.5*F[:,:,:,5];
mask = abs.(f0) .> 0 
fm = f0[mask]
fm = fm + randn(size(fm))/20
W = Diagonal(ones(N))
shat = -(W*H*A)\(W*f0)    % [9 1]. NB! May need to be rounded before applying settings on scanner.

f = f0 + H*A*shat;  % predicted fieldmap after applying shims

f0 = embed(f0, mask);
f = embed(f, mask);

% compare baseline and predicted fieldmaps
subplot(121)
im(f0); 
%title(sprintf('fieldmaps: original (left in each panel),\nand after 2nd order shimming (right)'));
title(sprintf('before shimming'));
h = colorbar; h.TickLabels{end} = 'Hz'; % h.Label.String = 'B0 field (Hz)';
subplot(122)
im(f,10*[-1 1]); 
title(sprintf('Residual (shim term %d)', shimNum));
h = colorbar; h.TickLabels{end} = 'Hz';
