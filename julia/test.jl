f(x,y,z) = [x y z]
@show f(1,2,3)
@show map(f, 1, 2, 3)
@show map(f, [1 2], [2 3], [3 4])
@show map(f, [1, 2], [2, 3], [3, 4])

using CoordinateTransformations: SphericalFromCartesian

	c2s = SphericalFromCartesian()

	function cart2sph(x,y,z)
		a = c2s([x,y,z]);
		(a.r, a.θ, π/2-a.ϕ)   # [radius azimuth polar-angle]
	end

#	x = [0
#	@show map(cart2sph(x,y,z), x, y, z)
