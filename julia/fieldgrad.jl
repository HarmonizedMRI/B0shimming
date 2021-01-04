using Interpolations

function fieldgrad(fov::Vector, f::Array)

	(nx,ny,nz) = size(f)

	itp = LinearInterpolation((1:nx,1:ny,1:nz), f)
	#itp = interpolate(f, BSpline(Quadratic(Line(OnGrid()))));
	
	g = [Interpolations.gradient(itp, x, y, z) for x in 1:nx, y in 1:ny, z in 1:nz]

	(dx, dy, dz) = fov ./ [nx,ny,nz]    # voxel size
	gx = map( g -> g[1]/dx, g);   # Hz/cm
	gy = map( g -> g[2]/dy, g);
	gz = map( g -> g[3]/dz, g);

	return (gx, gy, gz)
end

