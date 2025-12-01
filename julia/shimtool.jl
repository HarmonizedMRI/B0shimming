using MRIFieldmaps: b0map, b0model, b0init, b0scale, phasecontrast
using MIRTjim: jim, prompt; jim(:prompt, false)
using MIRTio: loadrds
using Statistics: mean
using Unitful: Hz, s
using Plots; default(markerstrokecolor=:auto, label="")
using ROMEO: unwrap
using NIfTI: niread, niwrite, NIVolume
using MAT, JLD2
using MIRT: embed!
using FFTW

#include("loadrds.jl")
include("Shimming/getSHbasis.jl")     # getSHbasis()
include("Shimming/getcalmatrix.jl")   # getcalmatrix()
include("Shimming/shimoptim.jl")      # shimoptim()
include("phant_mask.jl")              # phant_mask()
include("shim_mask.jl")               # shim_mask()
include("mask_THR.jl")                # mask_THR()
include("bs_mask.jl")                 # bs_mask()
include("fmap_pos.jl")                # fmap_pos()



#####

# Set Needed Parameters
echo = [2.22, 4.45];
ncoil = 32              # Number of coils
l2b = -3                # Regularization parameter
precon = :diag;         # Precondition
l = 2                   # 1 for linear shims, 2 for second order
TH = 0.4                # Default threshold for mask
fov = [24 24 24]        # Field of view in cm
nx = 60                 # Matrix Size
ny = 60                 # Matrix Size
nz = 60                 # Matrix Size
nz_dummy = 1            # Number of Dummy slices
mask_type = 2           # Mask type: 1. Threshold   2. Brain Extraction Tool  3. Brainstem  4. Weighted BS  5. Custom
fieldmatch = 0			# 1 To turn on image registration with previous data for field map matching
calibration_file = "shimcal1025_outGE.mat"
data_file = "path/to/file"

#####

clim = (-60, 60).*Hz
necho = length(echo)    # Number of Echos
echo = echo * 1f-3      # Echo time in seconds
mask_flag = 0           # Initalization of mask flag
nFramesPerZ = necho*ny  # Y dimension size (echos x ny)

# Open File 
fid = open(data_file)
dat = Array{Complex{Int16}}(undef, nx, ncoil, necho*ny, nz)
for z in nz_dummy .+ (1:nz)
    for y in 1:necho*ny
        frame = (z-1)*nFramesPerZ + y
        dat[:,:,y,z-nz_dummy] = loadrds(fid, frame, nx, ncoil);
    end
end
close(fid)

dat = permutedims(dat, [1, 3, 4, 2])     # Complex Int16, [nx nTE*ny nz ncoil]
dat = convert(Array{Complex{Float32}}, dat); 


##### Create Mask of ROI #####

# Create Magnitude Image
d1 = dat[:,1:2:end,:,:]
d2 = dat[:,2:2:end,:,:]
img1r = zeros(Complex,size(d1))
img2r = zeros(Complex,size(d2))

for ic = 1:ncoil
	img1r[:,:,:,ic] = fftshift(ifft(fftshift(d1[:,:,:,ic])));   # First echo multicoil image
	img2r[:,:,:,ic] = fftshift(ifft(fftshift(d2[:,:,:,ic])));   # Second echo multicoil image
end

img1 = sqrt.(sum(abs.(img1r),dims = 4));                        # First echo combined coil Image
mask=[]
mask_w=[]
maskbs=[]
maskbrain=[]
mask, mask_w,maskbs,maskbrain = shim_mask(img1, mask_type)    # Create Mask: shim_mask(magnitude image, mask type)


##### Create B0 Map #####

# Combine and scale complex data
images = cat(img1r, img2r; dims = 5)
img_sc = sum(images; dims=4)
(img_scale, scale) = b0scale(img_sc, echo)
images = images/scale;

# B0 Iniitalization and Unwrap
init1 = b0init(images,echo;)                                     # Returns Inital Map in Hz
init = init1.*mask                                               # Mask inital fieldmap
init_phase = (echo[2]-echo[1])*init*(2*pi)                       # Convert to Radians
uw_init1 = unwrap(init_phase[:,:,:,1]; mag = img1, mask = mask)  # Unwrap Returns Radians
uw_init = uw_init1 /(2*pi*(echo[2]-echo[1]))

# Regualrized B0 map
fmap = b0map(images,echo; smap = nothing, l2b, precon, finit = uw_init, mask = mask)[1]
pp = jim(fmap; title = "Field Map")
display(pp)


##### Shim calculations #####

loss = (s, HA, f0) -> norm(HA*s + f0, 2)^2 / length(f0)
ftol_rel = 0.5e-5

# Load Calibration Data
matf = matread(calibration_file);
F = matf["F"]               # [nx ny nz 8] for 2nd order shim systems
S = matf["S"]               # [8 8] matrix with shim amplitudes (typically diagonal)
mask_c = matf["mask_c"];    
mask_c = BitArray(mask_c)   # Calibration mask
fov = matf["FOV_c"];        # Calibration field of view

# Adjustment for first order (if needed)
if l < 2
	F = F[:,:,:,1:3]
	ss = diag(S)
	S = Diagonal(ss[1:3])
end

# Adjust Calibrtion data to fit shim order
inds = [3, 1, 2, 4, 6, 8, 7, 5]
Fr = copy(F)
for ii = 1:size(F,4)
	Fr[:,:,:,ii] = F[:,:,:,inds[ii]]
end
F = Fr
N = sum(vec(mask_c))        # Number of pixels in calibration mask
(nx,ny,nz,nShim) = size(F)
(x,y,z) = LinRange.(1, -1, [nx,ny,nz]) .* vec(fov)/2  # Spatial positions of pixels

# mask F and reshape to [N 8]
Fm = zeros(N, nShim)
for ii = 1:nShim
	f1 = F[:,:,:,ii]
	Fm[:,ii] = f1[mask_c]
end

# Get spherical harmonic basis of degree l for calibration data
H = getSHbasis(x, y, z; L=l) # [nx ny nz numSH(l)]
H = reshape(H, :, size(H,4))
H = H[vec(mask_c), :] # size is [N sum(2*(0:l) .+ 1)]

# Get calibration matrix A
A = getcalmatrix(Fm, H, diag(S))
f0m = fmap[mask]
N = sum(vec(mask))

# Get spherical harmonic basis of degree l for fieldmap data
H = getSHbasis(x, y, z; L=l)
H = reshape(H, :, size(H,4))
H = H[vec(mask), :]

# Determine weight (W)
if mask_w ==[]
    W = Diagonal(ones(N,))
else
    W = Diagonal(mask_w[mask])
end


if fieldmatch ==1
	fmapref = fmap_registration(img1)
end

# Solve for shim changes
s0 = -(W*H*A)\(W*f0m)    # Unconstrained least-squares solution
shat = s0


##### Output Shim Changes #####

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

println(" ")# Loss (objective) function for optimization.
# The field map f is modeled as f = H*A*s + f0, where
#   s = shim amplitudes (vector),
#   H = spherical harmonic basis functions
#   A = matrix containing shim coil expansion coefficients for basis in H
#   f0 = baseline field map at mask locations (vector)
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

# Predicted fieldmap after applying shims
fp = zeros(size(fmap))
fpm = H*A*shat + f0m
embed!(fp, fpm, mask)

# display predicted fieldmap
p = jim(log.(abs.(A[:,:]')); color=:RdBu)
p = jim(fp; clim=(-200,200), color=:RdBu)
p = jim(cat(fmap[:,:,:],fp[:,:,:];dims=1); ncol=6, clim=(-200,200), color=:RdBu)
display(p)

# predicted fieldmap after applying shims, no mask; Not necessary but useful for RMS comparisons
fpnm = zeros(size(fmap));
H1 = getSHbasis(x, y, z; L=l);
H1 = reshape(H1, :, size(H1,4));
fpnmv = H1*A*shat; 
embed!(fpnm, fpnmv, BitArray(ones(nx,ny,nz)));
fpnm = fpnm + fmap;                                 # Predicted fieldmap no mask 

rm_bet = (`rm /home/jpothoof/Julia_Files/bet_term_mask.nii`)
run(rm_bet)
println("Done")

# Output shims in text file
outfile = "shim_update.txt"
f = open(outfile, "w")
ge_string =  "$(shat_ge[1])" * " " * "$(shat_ge[3])" * " " * "$(shat_ge[4])" * " " * "$(shat_ge[2])" * " " * "$(shat_ge[5])" * " " * "$(shat_ge[6])" * " " * "$(shat_ge[7])" * " " * "$(shat_ge[8])" * " " * "$(shat_ge[9])" * " " * ""
println(f, ge_string)
close(f)
