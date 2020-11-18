
using CoordinateTransformations: SphericalFromCartesian
using SphericalHarmonics: computeYlm
using MIRT: jim, ndgrid

"""
	function getSHbasis(x,y,z,l)

Get spherical harmonic basis up to degree 'l' evaluated at [x y z]

Inputs:
  x/y/z    length-N vector     x/y/z coordinates at which to evaluate (discretize) the basis (cm)
  l        int                 SH degree. Default: 2.

Output:
  H       [N ...]    SH basis
"""
function getSHbasis(
	x::Vector{<:Real},
	y::Vector{<:Real},
	z::Vector{<:Real},
	L::Int = 2
)

	# Function to evaluate spherical harmonic of (degree,order) = (l,m) at position (x,y,z)
	c2s = SphericalFromCartesian()
	function sh(x,y,z,l,m)
		a = c2s([x,y,z]);
		(r,ϕ,θ) = (a.r, a.θ, π/2-a.ϕ)    # (radius, azimuth, colatitude)
		ylm = computeYlm(θ, ϕ, lmax = l)
		r^l * ylm[(l,m)]
	end

	# construct basis matrix H
	H = zeros(size(x,1), sum(2*(0:L) .+ 1));
	ic = 1;
	for l = 0:L
		for m = 0:l
			f = map((x,y,z) -> sh(x,y,z,l,m), x, y, z)
			H[:,ic] = real(f);
			ic = ic+1;
			if m > 0
				H[:,ic] = imag(f);
				ic = ic+1;
			end
 		end
	end

	# normalize DC column
	H[:,1] .= 1.

	H
end

function meshgrid(x, y, z)
	X = [i for i in x, j in y, k in z]
	Y = [j for i in x, j in y, k in z]
	Z = [k for i in x, j in y, k in z]
	return X, Y, Z
end

# test function
function getSHbasis(str::String)

	# define spatial locations
	r = range(-10,10,length=40)
	(x, y, z) = ndgrid(r, r, range(-10,10,length=20))

	# get spherical harmonic basis
	H = getSHbasis(x[:], y[:], z[:], 2)

	# reshape so each basis vector can be view with jim
	nb = size(H,2)
	(nx, ny, nz) = size(x)
	Hc = zeros(nx, ny, nz, nb)
	for ii = 1:nb
		Hc[:,:,:,ii] = reshape(H[:,ii], nx, ny, nz)
	end
	display(jim(Hc[:,:,:,2]))

	# save to .mat file for comparison with getSHbasis.m
	# matwrite("tmp.mat", Dict([("Hc", Hc), ("x", x), ("y", y), ("z", z)]))

	Hc    # return value
end
