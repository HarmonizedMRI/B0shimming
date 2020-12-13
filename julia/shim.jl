using MAT, JLD2
using MIRT: embed!, jim, ndgrid
#using LinearAlgebra
#using SparseArrays

include("getSHbasis.jl")
include("getcalmatrix.jl")
include("shimoptim.jl")


############################################################################################
## EDIT this section

# Shim calibration data
calFile = "CalibrationDataSiemensMGH12Dec2020.jld2"
shimlims = (100, 1000, 4000)
calFile = "CalibrationDataUM10Dec2020.jld2"
shimlims = (100, 4000, 12000)   # (max linear shim current, max hos shim current, max total hos shim current)

# baseline field map (and fov, mask). See mat2jld2.jl.
f0File = "f0_jar.jld2"   
f0File = "f0_redhead_localmask.jld2"   
f0File = "f0.jld2"   

# order of spherical harmonic basis
l = 2

# Loss (objective) function for optimization.
# The field map model is f = HA*s + f0, where
# s = shim amplitudes (vector), 
# HA = a matrix with calibrated basis functions, 
# f0 = baseline field map

loss = (s, HA, f0) -> norm(HA*s + f0)^2

############################################################################################

if l < 2
	F = F[:,:,:,1:3]
	s = diag(S)
	S = Diagonal(s[1:3])
end


############################################################################################
# Load calibration data and calculate calibration matrix A
############################################################################################

# F = [nx ny nz 8] (Hz), in order 'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'
# S = [8 8], shim amplitudes used to obtain F (hardware units)
@load "$calFile" F S fov mask

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
	range(-fov[1]/2, fov[1]/2, length=nx), 
	range(-fov[2]/2, fov[2]/2, length=ny), 
	range(-fov[3]/2, fov[3]/2, length=nz) 
	)

# mask F and reshape to [N 8] 
Fm = zeros(N, nShim)
for ii = 1:nShim
	f1 = F[:,:,:,ii]
	Fm[:,ii] = f1[mask]
end

# Get spherical harmonic basis of degree l
H = getSHbasis(x[mask], y[mask], z[mask], l)   # size is [N sum(2*(0:l) .+ 1)]

# Get calibration matrix A
A = getcalmatrix(Fm, H, diag(S))



############################################################################################
# Now that A is determined we can use it to optimize shims for a given acquired fieldmap f0.
############################################################################################

@load "$f0File" f0 fov mask

(nx,ny,nz) = size(f0)

mask = BitArray(mask)

f0m = f0[mask]

N = sum(mask[:])

(x,y,z) = ndgrid(
	range(-fov[1]/2, fov[1]/2, length=nx), 
	range(-fov[2]/2, fov[2]/2, length=ny), 
	range(-fov[3]/2, fov[3]/2, length=nz) 
	)

H = getSHbasis(x[mask], y[mask], z[mask], l)  
W = Diagonal(ones(N,))   # optional spatial weighting

s0 = -(W*H*A)\(W*f0m)    # Unconstrained least-squares solution (for comparison)
#@show Int.(round.(s0))

shat = shimoptim(W*H*A, W*f0m, shimlims; loss=loss) #, s0=[s0[1]; zeros(8,)])
@show Int.(round.(shat))

@show loss(zeros(9,), W*H*A, W*f0m)   # loss before optimization
@show loss(s0, W*H*A, W*f0m)          # loss after unconstrained optimization
@show loss(shat, W*H*A, W*f0m)        # loss after constrained optimization

shat_ge = Int.(round.(shat))
shat_siemens = round.(shat; digits=1)

println("\nRecommended shim changes:") 
println(string(
	"\tGE: cf, x, y, z = ", 
	shat_ge[1], ", ",
	shat_ge[3], ", ", 
	shat_ge[4], ", ", 
	shat_ge[2], " (set in Manual Prescan)")) 
println(string(
	"\tGE: setNavShimCurrent",
	" z2 ", shat_ge[5],
	" zx ", shat_ge[6], 
	" zy ", shat_ge[7], 
	" x2y2 ", shat_ge[8], 
	" xy ", shat_ge[9]
	) 
)
println(string(
	"\tSiemens: adjvalidate -shim -set -mp -delta ",
	shat_siemens[3], " ",
	shat_siemens[4], " ",
	shat_siemens[2], " ",
	shat_siemens[5], " ",
	shat_siemens[6], " ",
	shat_siemens[7], " ",
	shat_siemens[8], " ",
	shat_siemens[9]
	)
)

# predicted fieldmap after applying shims
fp = zeros(size(f0))
fpm = H*A*shat + f0m;      
embed!(fp, fpm, mask)   

# display predicted fieldmap
p = jim(log.(abs.(A[:,:]')); color=:jet)
p = jim(fp; clim=(-200,200), color=:jet)
p = jim(cat(f0,fp;dims=1); clim=(-200,200), color=:jet)    # compare before and after shimming
display(p)


# write to .mat file for viewing
matwrite("result.mat", Dict(
	"mask" => mask,
	"f0" => f0,
	"fp" => fp
))
