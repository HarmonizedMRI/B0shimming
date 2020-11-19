using JLD2
using MIRT: embed!, jim, ndgrid
using Random
using LinearAlgebra
using SparseArrays

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
A = getcalmatrix(Fm, H, S)

# Now that A it determined we can use it to optimize shims for a given acquired 'baseline' fieldmap.
# Example: Synthesize an example fieldmap 'f0' and optimize shims (minimize RMS residual) for that fieldmap.
f0 = F[:,:,:,2] + sqrt.(abs.(F[:,:,:,5]));  
f0 = sum(F[:,:,:,[2,3,8]], dims=4)[:,:,:,1]
mask = abs.(f0) .> 0                # note the dots
fm = f0[mask]
N = sum(mask[:])
fm = fm + 0*randn(size(fm))/20        # add some noise
W = sparse(collect(1:N), collect(1:N), ones(N))
shat = -(W*H*A)\(W*fm)    # Vector of length 9. NB! May need to be rounded before applying settings on scanner.

fpm = fm + H*A*shat;      # predicted fieldmap after applying shims

# display predicted fieldmap
fp = zeros(size(f0))
embed!(fp, fpm, mask)   
p = jim(cat(f0,fp;dims=1))    # compare before and after shimming
p = jim(fp; clim=(-20,20))
display(p)


