using JLD2
using MIRT: embed!, jim, ndgrid
using Random
using LinearAlgebra
using SparseArrays

include("getSHbasis.jl")
include("getcalmatrix.jl")
#include("ls_adam.jl")
include("shimoptim.jl")

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
mask = fm .> 200    # note the '.' (broadcasting)
mask[1:2:end, 1:2:end, 1:2:end] .= false
N = sum(mask[:])

# mask and reshape to [N 8] 
Fm = zeros(N, nShim)    # m for 'masked'
for ii = 1:nShim
	f1 = F[:,:,:,ii]
	Fm[:,ii] = f1[mask]
end

# Get spherical harmonic basis of degree l
l = 4
H = getSHbasis(x[mask], y[mask], z[mask], l)   # size is [N sum(2*(0:l) .+ 1)]

# Get calibration matrix A
A = getcalmatrix(Fm, H, S)

# Now that A is determined we can use it to optimize shims for a given acquired fieldmap.
# Example: Synthesize an example fieldmap 'f0' and optimize shims (minimize RMS residual) for that fieldmap.
f0 = sum(F[:,:,:,[2,6,7,8]], dims=4)[:,:,:,1]
f0 = F[:,:,:,2] + sqrt.(abs.(F[:,:,:,5])) + 3*F[:,:,:,6]
f0 = 2*F[:,:,:,2] + 10*sqrt.(abs.(F[:,:,:,5])) + 1*F[:,:,:,6]
f0 = F[:,:,:,2] + 2*F[:,:,:,6]
mask = abs.(f0) .> 0                # note the dots
mask[1:2:end, 1:2:end, 1:2:end] .= false
f0m = f0[mask]
N = sum(mask[:])
f0m = f0m + 0*randn(size(f0m))        # add some noise
H = getSHbasis(x[mask], y[mask], z[mask], l)  
W = sparse(collect(1:N), collect(1:N), ones(N))
W = Diagonal(ones(N,))

# Initial guess (unconstrained LS solution)
s0 = -(W*H*A)\(W*f0m)    # Vector of length 9. NB! May need to be rounded before applying settings on scanner.

# shim limits 
shimlims = (100, 2000, 12000)

# solve using Flux
# initialize with LS solution
(lin_max, hos_max, hos_sum_max) = shimlims
#shat[5:9] .= tan.(pi/2 * shat[5:9]/hos_max)   # initial guess
#@time (shat, out) = ls_adam(W*H*A, W*f0m; s0=s0, opt=opt, niter=niter)

# solve using NLopt
@time shat = shimoptim(W*H*A, W*f0m, shimlims) # ; s0=s0)

fpm = f0m + H*A*shat;      # predicted fieldmap after applying shims

# display predicted fieldmap
fp = zeros(size(f0))
embed!(fp, fpm, mask)   
p = jim(cat(f0,fp;dims=1))    # compare before and after shimming
p = jim(fp; clim=(-40,40))
#display(p)


