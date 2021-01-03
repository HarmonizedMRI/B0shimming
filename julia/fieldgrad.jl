

function fieldgrad(fov::Vector, f::Array)

	(nx,ny,nz) = size(f0)

	(dx, dy, dz) = fov ./ [nx,ny,nz]    # voxel size

	gx = zeros(size(f))
	gtmp = diff(f, dims=1) / dx
	gx[1:(end-1),:,:] = gtmp

	gy = zeros(size(f))
	gtmp = diff(f, dims=2) / dy
	gy[:,1:(end-1),:] = gtmp

	gz = zeros(size(f))
	gtmp = diff(f, dims=3) / dz
	gz[:,:,1:(end-1)] = gtmp

	return (gx, gy, gz)
end
