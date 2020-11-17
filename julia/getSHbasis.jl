
using CoordinateTransformations: SphericalFromCartesian
using Jacobi: legendre

"""
	function getSHbasis(x,y,z,l)

Evaluate spherical harmonic basis functions up to degree 'l' at locations [x y z]

Inputs:
  x/y/z    length-N vector     x/y/z coordinates at which to evaluate (discretize) the basis (cm)
  l        int                 SH degree. Default: 4.

Output:
  H       [N ...]    SH basis, including dc (B0 offset) term
"""
function getSHbasis(
	x::Vector{<:Real},
	y::Vector{<:Real},
	z::Vector{<:Real},
	L::Int = 4
)

	# Function evaluates spherical harmonic of (degree,order) = (l,m) at position (x,y,z)
	c2s = SphericalFromCartesian()
	function sh(x,y,z)              # l,m defined outside function
		a = c2s([x,y,z]);
		r = a.r                      # radius
		ϕ = a.θ                      # azimuth
		θ = π/2-a.ϕ                  # polar angle
		r^l * exp(im*m*ϕ) * legendre(cos(θ),l)
	end

	# construct basis matrix H
	H = zeros(size(x,1), sum(2*(0:L)+1));
	ic = 1;
	for l = 0:L
		for m = -l:l
			@show f = map(sh, x, y, z)
			#H(:,ic) = real(f);
			#ic = ic+1;
			#if m > 0
			#	H(:,ic) = imag(f);
			#	ic = ic+1;
 			#end
 		end
	end

end

