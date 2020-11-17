using CoordinateTransformations: SphericalFromCartesian

"""
	cart2sph(x,y,z)

Concvert vectors x/y/z of cartesian coordinates to
vectors r/ϕ/θ of spherical coordinates (radius, azimuth, polar angle)

Spherical coordinate system as defined in 
https://en.wikipedia.org/wiki/Spherical_harmonics
and Romeo and Hoult MRM 1984


"""
function cart2sph(
	x::Vector{<:Real},
	y::Vector{<:Real},
	z::Vector{<:Real},
)

	c2s = SphericalFromCartesian()

	function getr(x,y,z)
		a = c2s([x,y,z]);
		#(a.r, a.θ, π/2-a.ϕ)   # [radius azimuth polar-angle]
		a.r
	end

	function getϕ(x,y,z)
		a = c2s([x,y,z]);
		a.θ                    # azimuth
	end

	function getθ(x,y,z)
		a = c2s([x,y,z]);
		π/2-a.ϕ                # polar angle
	end

	r = map((x,y,z) -> getr(x,y,z), x, y, z)
	ϕ = map((x,y,z) -> getϕ(x,y,z), x, y, z)
	θ = map((x,y,z) -> getθ(x,y,z), x, y, z)

	(r,ϕ,θ)

end
