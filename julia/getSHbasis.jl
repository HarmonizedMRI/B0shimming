
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
	a = cart2sph(x, y, z)
	r = a[1]

	N = length(x)
	p = Vector{Vector{Float64}}(Float64, N)
	for ii = 1:N
		p[ii] = [x(ii], y(jj), z(jj)]
	end

	a = map((r,g) -> signalloss(r,g,Δ,A,Δs,te), r, g)
	
	%r = SMatrix{N,3}([x y z])

	c2s = SphericalFromCartesian()

	function cart2sph(x,y,z)
		a = c2s(p);      # a = (radius, azimuth, elevation)
		r = a.r          # radius
		ϕ = a.θ          # azimuth
		θ = π/2 - a.ϕ    # polar angle 
		[r, ϕ, θ]
	end

	sc = map( (x,y,z) -> cart2sph([x,y,z])
		
		
	end
	f(x,y,z)

	# Spherical coordinate system as defined in https://en.wikipedia.org/wiki/Spherical_harmonics
	# and Romeo and Hoult MRM 1984

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

# test function
function getSHbasis(str::String)
	x = [1,0,0]
	y = [0,1.,0]
	z = [0,0,1.]
)

end
