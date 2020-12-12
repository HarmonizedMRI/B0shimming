using MAT, JLD2
using MIRT: embed!, jim, ndgrid
#using LinearAlgebra
#using SparseArrays

include("getSHbasis.jl")
include("getcalmatrix.jl")
include("shimoptim.jl")

############################################################################################
# Load calibration data and calculate calibration matrix A
############################################################################################

# F = [nx ny nz 8] (Hz), in order 'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'
# S = [8 8], shim amplitudes used to obtain F (hardware units)
@load "CalibrationDataUM10Dec2020.jld2" F S fov mask
mask = BitArray(mask)

# 0th-2nd order terms in getSHbasis.jl are in order [dc z x y z2 zx zy x2y2 xy],
# so reorder F to match that.
# No need to reorder S
inds = [3, 1, 2, 4, 6, 8, 7, 5] 
Fr = copy(F)
for ii = 1:size(F,4)
	Fr[:,:,:,ii] = F[:,:,:,inds[ii]]
end
F = Fr

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
l = 2
H = getSHbasis(x[mask], y[mask], z[mask], l)   # size is [N sum(2*(0:l) .+ 1)]

# Get calibration matrix A
#S = collect(Diagonal(ones(8,)))
A = getcalmatrix(Fm, H, S)
#@show size(A)



############################################################################################
# Now that A is determined we can use it to optimize shims for a given acquired fieldmap f0.
############################################################################################

@load "f0.jld2" f0 fov mask

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

# define loss 
loss = (s, HA, f0) -> norm(HA*s + f0)^2

s0 = -(W*H*A)\(W*f0m)    # Unconstrained least-squares solution
@show Int.(round.(s0))

shat = shimoptim(W*H*A, W*f0m, shimlims; loss=loss) #, s0=[s0[1]; zeros(8,)])
@show Int.(round.(shat))
#shat = s0

@show loss(zeros(9,), W*H*A, W*f0m)
@show loss(s0, W*H*A, W*f0m)
@show loss(shat, W*H*A, W*f0m)

shat = Int.(round.(shat))

println("\nRecommended shim changes:") 
println(string(
	"\tcf, z, x, y = ", 
	-shat[1], ", ",    # whether the minus sign is needed here is likely scanner-specific
	shat[2],  ", ", 
	shat[3],  ", ", 
	shat[4],  " (set in Manual Prescan)")) 
println(string(
	"\tsetNavShimCurrent z2 ", shat[5],
	" zx ", shat[6], 
	" zy ", shat[7], 
	" x2y2 ", shat[8], 
	" xy ", shat[9]
	) 
)

# predicted fieldmap after applying shims
fp = zeros(size(f0))
fpm = H*A*shat + f0m;      
embed!(fp, fpm, mask)   

# display predicted fieldmap
p = jim(log.(abs.(A[:,:]')); color=:jet)
p = jim(cat(f0,fp;dims=1); clim=(-80,80), color=:jet)    # compare before and after shimming
p = jim(fp; clim=(-50,50), color=:jet)
display(p)


# write to .mat file for viewing
matwrite("result.mat", Dict(
	"fp" => fp
))
