using JLD2
using MIRT: embed!, jim, ndgrid
using Random
using LinearAlgebra
using SparseArrays

include("getSHbasis.jl")
include("getcalmatrix.jl")
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

# mask F and reshape to [N 8] 
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
f0 = 2*F[:,:,:,2] + 10*sqrt.(abs.(F[:,:,:,5])) + 0.3*F[:,:,:,6]
f0 = 2*F[:,:,:,2] + 1*F[:,:,:,6]
mask = abs.(f0) .> 0                # note the dots
mask[1:2:end, 1:2:end, 1:2:end] .= false
f0m = f0[mask]
N = sum(mask[:])
#f0m = f0m + 0*randn(size(f0m))        # add some noise

H = getSHbasis(x[mask], y[mask], z[mask], l)  
W = sparse(collect(1:N), collect(1:N), ones(N))
W = Diagonal(ones(N,))

# Unconstrained LS solution
# s0 = -(W*H*A)\(W*f0m)

# shim limits 
shimlims = (100, 4000, 12000)   # (max linear shim current, max hos shim current, max total hos shim current)

# define loss and solve for shims 
loss = (s, HA, f0) -> norm(HA*s + f0)^2
@time shat = shimoptim(W*H*A, W*f0m, shimlims; loss=loss) #;  s0=s0)

# return integer values
s = Int.(round.(shat))
println("Optimized shims: $s")

# predicted fieldmap after applying shims
fpm = H*A*s + f0m;      

# display predicted fieldmap
fp = zeros(size(f0))
embed!(fp, fpm, mask)   
p = jim(fp; clim=(-40,40))
p = jim(cat(f0,fp;dims=1))    # compare before and after shimming
display(p)


