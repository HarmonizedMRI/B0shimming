# Julia script for calculating optimal shim settings
#
# Run this function from the HarmonizedMRI/B0shimming/julia directory
#
# This script requires the following files to exist in the current path:
#   shimcal.mat      F S mask FOV
#   f0.mat           f0 FOV. Defines N and FOV used to create shimvol.mat.
#   shimvol.mat      mask. Shim ROI mask on grid defined by N and FOV,
#                    excluding points outside object mask (defined in f0.mat)
#
# In addition, to use SlicePlanner to select an ROI, you'll need:
#   Localizer.h5     For displaying object in SlicePlanner GUI

using MAT, JLD2
using MIRT: embed!
using MIRTjim: jim

include("getSHbasis.jl")     # getSHbasis()
include("getcalmatrix.jl")   # getcalmatrix()
include("shimoptim.jl")      # shimoptim()

# shim system current limits
shimlims = (5000*ones(3,), 4000*ones(5,), 12000)   # (lin max, hos max, sum hos max). For GE 3T scanner at U of Michigan.

# Order of spherical harmonic basis.
# For linear shimming, set l = 1.
l = 2

# Loss (objective) function for optimization.
# The field map f is modeled as f = H*A*s + f0, where
#   s = shim amplitudes (vector),
#   H = spherical harmonic basis functions
#   A = matrix containing shim coil expansion coefficients for basis in H
#   f0 = baseline field map at mask locations (vector)
loss = (s, HA, f0) -> norm(HA*s + f0, 2)^2 / length(f0)

ftol_rel = 0.5e-5


############################################################################################
# Calculate calibration matrix A
############################################################################################

# For 2nd order shims:
# F = [nx ny nz 8] (Hz), in order 'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'
# S = [8 8], shim amplitudes used to obtain F (hardware units)
matf = matread("shimcal.mat");
F = matf["F"]   # [nx ny nz 8] for 2nd order shim systems
S = matf["S"]   # [8 8] matrix with shim amplitudes (typically diagonal)
mask = matf["mask"];
fov = matf["FOV"];

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

N = sum(vec(mask))

(nx,ny,nz,nShim) = size(F)

(x,y,z) = LinRange.(1, -1, [nx,ny,nz]) .* vec(fov)/2

# mask F and reshape to [N 8]
Fm = zeros(N, nShim)
for ii = 1:nShim
	f1 = F[:,:,:,ii]
	Fm[:,ii] = f1[mask]
end

# Get spherical harmonic basis of degree l
H = getSHbasis(x, y, z; L=l) # [nx ny nz numSH(l)]
H = reshape(H, :, size(H,4))
H = H[vec(mask), :] # size is [N sum(2*(0:l) .+ 1)]

# Get calibration matrix A
A = getcalmatrix(Fm, H, diag(S))

# @save "A_l$l.jld2" A # uncomment this if you want to save it



############################################################################################
# Now that A is determined we can use it to optimize shims for a given acquired fieldmap f0.
############################################################################################

matf = matread("f0.mat");
f0 = matf["f0"]   # B0 field map we wish to shim
fov = matf["FOV"];

matf = matread("shimvol.mat")
mask = matf["mask"]

(nx,ny,nz) = size(f0)

mask = BitArray(mask)

f0m = f0[mask]

N = sum(vec(mask))

(x,y,z) = LinRange.(1, -1, [nx,ny,nz]) .* vec(fov)/2

H = getSHbasis(x, y, z; L=l)
H = reshape(H, :, size(H,4))
H = H[vec(mask), :]
W = Diagonal(ones(N,))   # optional spatial weighting

s0 = -(W*H*A)\(W*f0m)    # Unconstrained least-squares solution (for comparison)

# This is where it all happens.
# @time shat = shimoptim(W*H*A, W*f0m, shimlims; loss=loss, ftol_rel=ftol_rel)
shat = s0
@show Int.(round.(shat))

# Print and plot results

@show loss(0*s0, W*H*A, W*f0m)   # loss before optimization
@show loss(s0, W*H*A, W*f0m)     # loss after unconstrained optimization
@show loss(shat, W*H*A, W*f0m)   # loss after constrained optimization

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
	"\tset cf, x, y, z shims in Manual Prescan"))
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
fpm = H*A*shat + f0m
embed!(fp, fpm, mask)

# println()

# display predicted fieldmap
p = jim(log.(abs.(A[:,:]')); color=:jet)
p = jim(fp; clim=(-200,200), color=:jet)
iz = 16:45 # compare these slices before and after shimming:
p = jim(cat(f0[:,:,iz],fp[:,:,iz];dims=1); ncol=6, clim=(-200,200), color=:jet)
display(p)

# predicted fieldmap after applying shims, no mask
fpnm = zeros(size(f0));
H = getSHbasis(x, y, z; L=l);
H = reshape(H, :, size(H,4));
fpnmv = H*A*shat;
embed!(fpnm, fpnmv, BitArray(ones(nx,ny,nz)));
fpnm = fpnm + f0;

# Optional: write to .mat file for viewing
# requires MAT package
matwrite("result.mat", Dict(
	"mask" => mask,
	"f0" => f0,
	"fp" => fp,
	"fpnm" => fpnm
))
