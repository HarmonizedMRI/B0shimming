using SphericalHarmonics
using CoordinateTransformations: SphericalFromCartesian
using MIRT: ndgrid, jim
using ForwardDiff

c2s = SphericalFromCartesian()

# function evaluates spherical harmonic at one spatial location r = [x, y, z]
function c2sph(r, l, m)
	a = c2s(r)
	(rho, ϕ, θ) = (a.r, a.θ, π/2 - a.ϕ)    # (radius, azimuth, colatitude)
	Y = computeYlm(θ, ϕ; lmax=l)
	rho^l * Y[(l,m)]
end

"""
	H = getSHbasis(x::Vector, y::Vector, z::Vector, L::Int64)

Get spherical harmonic basis of order L evaluated at spatial locations x, y, z
Order of 0th-2nd order terms:  
	H[:,1]   cf (center frequency, Hz)  
	H[:,2]   z
	H[:,3]   x
	H[:,4]   y
	H[:,5]   z2
	H[:,6]   zx
	H[:,7]   zy
	H[:,8]   x2y2
	H[:,9]   xy
"""
function getSHbasis(
	x::Vector{<:Real}, 
	y::Vector{<:Real}, 
	z::Vector{<:Real}, 
	L::Int64
	)

	# vector of position vectors 
	r = map( (x,y,z) -> [x, y, z], x, y, z)

	# SH basis
	H = zeros(size(x,1), sum(2*(0:L) .+ 1))
	ic = 1
	for l = 0:L
		for m = -0:l
			# SH evaluated at r
			f = map( r -> c2sph(r, l, m), r)
			H[:,ic] = real(f)    

			ic = ic+1
			if m != 0
				H[:,ic] = imag(f)
				ic = ic+1
 			end
 		end
 	end

	H[:,1] .= 1.0    # center frequency offset

	H
end

function getSHbasisGrad(
	x::Vector{<:Real}, 
	y::Vector{<:Real}, 
	z::Vector{<:Real}, 
	L::Int64
	)

	# vector of position vectors 
	r = map( (x,y,z) -> [x, y, z], x, y, z)

	# SH basis
	dHx = zeros(size(x,1), sum(2*(0:L) .+ 1))
	dHy = zeros(size(x,1), sum(2*(0:L) .+ 1))
	dHz = zeros(size(x,1), sum(2*(0:L) .+ 1))
	ic = 1
	for l = 0:L
		for m = -0:l
			sh1 =  r -> real(c2sph(r, l, m))
			g = r -> ForwardDiff.gradient(sh1, r)
			df = map( r -> g(r), r)
			dHx[:,ic] = map( df -> df[1], df)
			dHy[:,ic] = map( df -> df[2], df)
			dHz[:,ic] = map( df -> df[2], df)

			ic = ic+1
			if m != 0
				sh1 =  r -> imag(c2sph(r, l, m))
				g = r -> ForwardDiff.gradient(sh1, r)
				df = map( r -> g(r), r)
				dHx[:,ic] = map( df -> df[1], df)
				dHy[:,ic] = map( df -> df[2], df)
				dHz[:,ic] = map( df -> df[2], df)

				ic = ic+1
 			end
 		end
 	end

	(dHx, dHy, dHz)
end

# test function
function getSHbasis(str::String)

	(nx,ny,nz) = (40,40,20)
	rx = range(-10,10,length=nx)
	ry = range(-10,10,length=ny)
	rz = range(-10,10,length=nz)
	(x,y,z) = ndgrid(rx,ry,rz)

	l = 2
	H = getSHbasis(x[:], y[:], z[:], l)
	(dHx, dHy, dHz) = getSHbasisGrad(x[:], y[:], z[:], l)

	nb = size(H,2)
	H = reshape(H, nx, ny, nz, nb)
	dHx = reshape(dHx, nx, ny, nz, nb)
	dHy = reshape(dHy, nx, ny, nz, nb)
	dHz = reshape(dHz, nx, ny, nz, nb)

	# jim(H[:,:,:,7], color=:jet)   # compare with >> evalspharm("test")

	(H, dHx, dHy, dHz)
end

