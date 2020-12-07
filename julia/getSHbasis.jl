using SphericalHarmonics
using CoordinateTransformations: SphericalFromCartesian
using MIRT: ndgrid

function getSHbasis(
	x::Vector{<:Real}, 
	y::Vector{<:Real}, 
	z::Vector{<:Real}, 
	L::Int64
	)

	# function evaluates spherical harmonic at one spatial location (x,y,z)
	c2s = SphericalFromCartesian()
	function sh(x, y, z, l, m)
		a = c2s([x,y,z])
		(r, ϕ, θ) = (a.r, a.θ, π/2 - a.ϕ)    # (radius, aziumuth, colatitude)
		Y = computeYlm(θ, ϕ; lmax=l)
		f = r^l * Y[(l,m)]
		return f
	end

	H = zeros(size(x,1), sum(2*(0:L) .+ 1))
	ic = 1
	for l = 0:L
		for m = 0:l
			f = map( (x,y,z) -> sh(x, y, z, l, m), x, y, z)
			H[:,ic] = real(f)
			ic = ic+1
			if m > 0
				H[:,ic] = imag(f)
				ic = ic+1
 			end
 		end
 	end

	H
end

function getSHbasis(str::String)

	(nx,ny,nz) = (40,40,20)
	rx = range(-10,10,length=nx)
	ry = range(-10,10,length=ny)
	rz = range(-10,10,length=nz)
	(x,y,z) = ndgrid(rx,ry,rz)

	l = 2
	H = getSHbasis(x[:], y[:], z[:], l)

	H = reshape(H, nx, ny, nz, size(H,2))
end
