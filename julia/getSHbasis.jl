using SphericalHarmonics
using CoordinateTransformations: SphericalFromCartesian

function getSHbasis(
	x::Vector{<:Real}, 
	y::Vector{<:Real}, 
	z::Vector{<:Real}, 
	L::Int64
	)

	# Spherical coordinate system as defined in https://en.wikipedia.org/wiki/Spherical_harmonics
	# and Romeo and Hoult MRM 1984
	c2s = SphericalFromCartesian()
	function sh(x, y, z, l, m)
		a = c2s([x,y,z])
		(r, ϕ, θ) = (a.r, a.θ, π/2 - a.ϕ)    # (radius, aziumuth, colatitude)
		Y = computeYlm(θ, ϕ; lmax=l)
		Ylm[(l,m)]
	end

	# construct basis matrix H
	H = zeros(size(x,1), sum(2*(0:L) .+ 1))
	ic = 1
	for l = 0:L
		for m = 0:l1
			f = map( (x,y,z) -> sh(x, y, z, l, m), x, y, z)
			H(:,ic) = real(f)
			ic = ic+1
			if m > 0
				H(:,ic) = imag(f)
				ic = ic+1
 			end
 		end
 	end
end

