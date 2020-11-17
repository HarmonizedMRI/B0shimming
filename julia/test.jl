function f(x,y,z) 
	x^2, y^2, z^2
end

@show f(1,2,3)
@show map(f, 1, 2, 3)
@show map(f, [1 2], [2 3], [3 4])
@show a = map(f, [1, 2], [2, 3], [3, 4])

using CoordinateTransformations: SphericalFromCartesian

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

	x = [1., 0, 0]
	y = [0., 1, 0]
	z = [0., 0, 1]

	r = map((x,y,z) -> getr(x,y,z), x, y, z)
	ϕ = map((x,y,z) -> getr(x,y,z), x, y, z)
	r = map((x,y,z) -> getr(x,y,z), x, y, z)

#	x = [0
#	@show map(cart2sph(x,y,z), x, y, z)
