using MAT, JLD2
using MIRT: embed!, jim, ndgrid
#using LinearAlgebra
#using SparseArrays

include("getSHbasis.jl")
include("getcalmatrix.jl")
include("shimoptim.jl")

"""
Load calibration data and calculate calibration matrix A
"""

@load "CalibrationDataUM10Dec2020.jld2" F S fov mask

mask = BitArray(mask)
N = sum(mask[:])

(nx,ny,nz,nShim) = size(F)

(x,y,z) = ndgrid(
	range(-fov[1], fov[1], length=nx), 
	range(-fov[2], fov[2], length=ny), 
	range(-fov[3], fov[3], length=nz) 
	)

# mask F and reshape to [N 8] 
Fm = zeros(N, nShim)
for ii = 1:nShim
	f1 = F[:,:,:,ii]
	Fm[:,ii] = f1[mask]
end

# Get spherical harmonic basis of degree l
l = 4
@time H = getSHbasis(x[mask], y[mask], z[mask], l)   # size is [N sum(2*(0:l) .+ 1)]

# Get calibration matrix A
#S = collect(Diagonal(ones(8,)))
A = getcalmatrix(Fm, H, S)


"""
Now that A is determined we can use it to optimize shims for a given acquired fieldmap f0.
"""

@load "f0.jld2" f0 fov mask

# f0 = F[:,:,:,2] + F[:,:,:,8]
mask = BitArray(mask)

f0m = f0[mask]

N = sum(mask[:])

(x,y,z) = ndgrid(
	range(-fov[1], fov[1], length=nx), 
	range(-fov[2], fov[2], length=ny), 
	range(-fov[3], fov[3], length=nz) 
	)

H = getSHbasis(x[mask], y[mask], z[mask], l)  
W = Diagonal(ones(N,))   # optional spatial weighting

# shim limits 
shimlims = (100, 4000, 12000)   # (max linear shim current, max hos shim current, max total hos shim current)

# define loss and 
loss = (s, HA, f0) -> norm(HA*s + f0)^2

s0 = zeros(9,)
@show loss(s0, W*H*A, W*f0m)
@time shat = shimoptim(W*H*A, W*f0m, shimlims; loss=loss,  s0=s0)
@show loss(shat, W*H*A, W*f0m)

s = Int.(round.(shat))
println("Optimized shims: $s")

# predicted fieldmap after applying shims
fp = zeros(size(f0))
fpm = H*A*s + f0m;      
embed!(fp, fpm, mask)   

# display predicted fieldmap
p = jim(cat(f0,fp;dims=1); clim=(-50,50), color=:jet)    # compare before and after shimming
display(p)

# write to .mat file for viewing
matwrite("result.mat", Dict(
	"f0" => f0,
	"fp" => fp
))
