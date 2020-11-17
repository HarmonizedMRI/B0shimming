
using StaticArrays: SVector
using CoordinateTransformations: SphericalFromCartesian

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
	l::Int = 4
)

	# convert to spherical coordinates
	N = length(x)
	p = Vector{Vector{Float64}}(Float64, N)
	for ii = 1:N
		p[ii] = [x(ii], y(jj), z(jj)]
	end
	
	%r = SMatrix{N,3}([x y z])

	cart2sph = SphericalFromCartesian()

	# Spherical coordinate system as defined in https://en.wikipedia.org/wiki/Spherical_harmonics
	# and Romeo and Hoult MRM 1984
	a = cart2sph.(p);     # a = (radius, azimuth, elevation)
	r = a.r          # radius
	θ = π/2 - a.ϕ    # polar angle 
	ϕ = a.θ          # azimuth

	# construct basis matrix H
	H = zeros(size(x,1), sum(2*(0:l)+1));
	ic = 1;
	for l1 = 0:l
		lp = legendre(l1, cos(θ));   % [l1+1 N]
		for m = 0:l1
			f = r.^l1 .* exp(1i*m*ϕ) .* lp(m+1,:)';
			H(:,ic) = real(f);
			ic = ic+1;
			if m > 0
				H(:,ic) = imag(f);
				ic = ic+1;
 			end
 		end
	end

end

function getSHbasis(
	x::Real,
	y::Real,
	z::Real,
	l::Int = 4
)

end
