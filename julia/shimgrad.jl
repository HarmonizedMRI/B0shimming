using MAT, JLD2
using MIRT: embed!, jim, ndgrid
#using LinearAlgebra
#using SparseArrays

include("getSHbasis.jl")
include("getcalmatrix.jl")
include("shimoptim.jl")
include("fieldgrad.jl")


############################################################################################
## EDIT this section

# Shim calibration data
calFile = "CalibrationDataSiemensMGH12Dec2020.jld2"
calFile = "CalibrationDataUM10Dec2020.jld2"; 

# shim system current limits
shimlims = (100*ones(3,), 4000*ones(5,), 12000)   # (lin max, hos max, sum hos max)

# baseline field map, fov, and mask. See mat2jld2.jl.
f0File = "Psub1.jld2"   

# order of spherical harmonic basis
# for linear shimming, set l = 1
l = 2

# Loss (objective) function for optimization.
# The field map f is modeled as f = H*A*s + f0, where
#   s = shim amplitudes (vector), 
#   H = spherical harmonic basis functions
#   A = matrix containing shim coil expansion coefficients for basis in H
#   f0 = baseline field map at mask locations (vector)

function loss(s, dHxA, dHyA, dHzA, g0x, g0y, g0z) 
	# TODO: account for non-isotropic voxel size
	g = map( (gx,gy,gz) -> norm([gx,gy,gz],2)^2, dHxA*s + g0x, dHyA*s + g0y, dHzA*s + g0z)
	return norm(g, 2)^2
end

ftol_rel = 1e-5

############################################################################################



############################################################################################
# Load calibration data and calculate calibration matrix A
############################################################################################

# For 2nd order shims: 
# F = [nx ny nz 8] (Hz), in order 'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'
# S = [8 8], shim amplitudes used to obtain F (hardware units)
@load "$calFile" F S fov mask

mask = BitArray(mask)

if l < 2
	F = F[:,:,:,1:3]
	s = diag(S)
	S = Diagonal(s[1:3])
end

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

(g0x,g0y,g0z) = fieldgrad(fov, f0);

g0xm = g0x[mask]
g0ym = g0y[mask]
g0zm = g0z[mask]

N = sum(mask[:])

(x,y,z) = ndgrid(
	range(-fov[1]/2, fov[1]/2, length=nx), 
	range(-fov[2]/2, fov[2]/2, length=ny), 
	range(-fov[3]/2, fov[3]/2, length=nz) 
	)

@time (dHx, dHy, dHz) = getSHbasisGrad(x[mask], y[mask], z[mask], l)  

W = Diagonal(ones(N,))   # optional spatial weighting

@time shat = shimoptim(W*dHx*A, W*dHy*A, W*dHz*A, g0xm, g0ym, g0zm, shimlims, loss)
@show Int.(round.(shat))

@show loss(0*shat, W*dHx*A, W*dHy*A, W*dHz*A, g0xm, g0ym, g0zm)   # loss before optimization
@show loss(  shat, W*dHx*A, W*dHy*A, W*dHz*A, g0xm, g0ym, g0zm)   # loss after contrained optimization

shat_ge = Int.(round.(shat))
shat_siemens = round.(shat; digits=1)

println("\nRecommended shim changes:") 

println(string(
	"\tcf, x, y, z = ", 
	shat_ge[1], ", ",
	shat_ge[3], ", ", 
	shat_ge[4], ", ", 
	shat_ge[2])) 

if length(shat) > 4
	println(string("\t",
		"z2 ", shat_ge[5],
		" zx ", shat_ge[6], 
		" zy ", shat_ge[7], 
		" x2y2 ", shat_ge[8], 
		" xy ", shat_ge[9]))

println(" ")
println(string(
	"GE: ",
	"\tset cf, x, y, z shims in Manual Prescan"));
	println(string(
		"\tsetNavShimCurrent",
		" z2 ", shat_ge[5],
		" zx ", shat_ge[6], 
		" zy ", shat_ge[7], 
		" x2y2 ", shat_ge[8], 
		" xy ", shat_ge[9]))
	println(" ")
	println(string(
		"Siemens: adjvalidate -shim -set -mp -delta ",
		shat_siemens[3], " ",
		shat_siemens[4], " ",
		shat_siemens[2], " ",
		shat_siemens[5], " ",
		shat_siemens[6], " ",
		shat_siemens[7], " ",
		shat_siemens[8], " ",
		shat_siemens[9]))
else
	println(string(
		"\tSiemens: adjvalidate -shim -set -mp -delta ",
		shat_siemens[3], " ",
		shat_siemens[4], " ",
		shat_siemens[2]))
end

# predicted fieldmap after applying shims
fp = zeros(size(f0))
H = getSHbasis(x[mask], y[mask], z[mask], l)   # size is [N sum(2*(0:l) .+ 1)]
f0m = f0[mask]
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
