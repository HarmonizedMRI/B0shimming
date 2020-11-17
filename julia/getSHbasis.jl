
using CoordinateTransformations: SphericalFromCartesian
using SphericalHarmonics: computeYlm

"""
	function getSHbasis(x,y,z,l)

Evaluate spherical harmonic basis functions up to degree 'l' at locations [x y z]

Inputs:
  x/y/z    length-N vector     x/y/z coordinates at which to evaluate (discretize) the basis (cm)
  l        int                 SH degree. Default: 2.

Output:
  H       [N ...]    SH basis, including dc (B0 offset) term
"""
function getSHbasis(
	x::Vector{<:Real},
	y::Vector{<:Real},
	z::Vector{<:Real},
	L::Int = 2
)

	c2s = SphericalFromCartesian()
	function sh(x,y,z,l,m)
		# Evaluate spherical harmonic of (degree,order) = (l,m) at position (x,y,z)
		a = c2s([x,y,z]);
		(r,ϕ,θ) = (a.r, a.θ, π/2-a.ϕ)    # (radius, azimuth, colatitude)
		ylm = computeYlm(θ, ϕ, lmax = l)
		# plm = computePlmcostheta(θ,lmax = l)
		#r^l * exp(-im*m*ϕ) * plm[(l,m)]
		r^l * ylm[(l,m)]
	end

	# construct basis matrix H
	H = zeros(size(x,1), sum(2*(0:L) .+ 1));
	ic = 1;
	for l = 0:L
		for m = -l:l
			f = map((x,y,z) -> sh(x,y,z,l,m), x, y, z)
			H[:,ic] = real(f);
			ic = ic+1;
 		end
	end

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
	(x, y, z) = meshgrid(
	x = [(i-nx/2)/(nx/2) for i = 1:nx, j = 1:ny, z = 1:nz]

end
