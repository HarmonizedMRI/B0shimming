using SphericalHarmonics
using CoordinateTransformations: SphericalFromCartesian
using MIRT: ndgrid, jim
using ForwardDiff
using AutoGrad
using ReverseDiff

# function evaluates spherical harmonic at one spatial location r = [x, y, z]
c2s = SphericalFromCartesian()
function c2sph(r, l, m)
	a = c2s(r)
	(rho, ϕ, θ) = (a.r, a.θ, π/2 - a.ϕ)    # (radius, azimuth, colatitude)
	Y = computeYlm(θ, ϕ; lmax=l)
	rho^l * Y[(l,m)]
end

"""
	H = getSHbasis(x, y, z, L)

Get spherical harmonic basis up to order L evaluated at spatial locations x, y, z
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

	r = map( (x,y,z) -> [x, y, z], x, y, z)

	H = zeros(size(x,1), sum(2*(0:L) .+ 1))

	ic = 1
	for l = 0:L
		for m = -0:l
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

	return H
end

"""
	getSHbasisGrad(x, y, z, L)

Get gradient of spherical harmonic basis up to order L evaluated at spatial locations x, y, z
"""
function getSHbasisGrad(
	x::Vector{<:Real}, 
	y::Vector{<:Real}, 
	z::Vector{<:Real}, 
	L::Int64
	)

	r = map( (x,y,z) -> [x, y, z], x, y, z)

	# dH = Array{Vector}(undef, length(x), sum(2*(0:L) .+ 1))
	dHx = zeros(size(x,1), sum(2*(0:L) .+ 1))
	dHy = zeros(size(x,1), sum(2*(0:L) .+ 1))
	dHz = zeros(size(x,1), sum(2*(0:L) .+ 1))

	ic = 1
	for l = 0:L
		for m = -0:l
			sh1 =  r -> real(c2sph(r, l, m))
			df = map( r -> ForwardDiff.gradient(sh1, r), r)
			dHx[:,ic] = map( df -> df[1], df)
			dHy[:,ic] = map( df -> df[2], df)
			dHz[:,ic] = map( df -> df[3], df)
			ic = ic+1
			if m != 0
				sh1 =  r -> imag(c2sph(r, l, m))
				df = map( r -> ForwardDiff.gradient(sh1, r), r)
				dHx[:,ic] = map( df -> df[1], df)
				dHy[:,ic] = map( df -> df[2], df)
				dHz[:,ic] = map( df -> df[3], df)
				ic = ic+1
 			end
 		end
 	end

	dHx[isnan.(dHx)] .= 0;
	dHy[isnan.(dHy)] .= 0;
	dHz[isnan.(dHz)] .= 0;

	return (dHx, dHy, dHz)
end


"""
	getSHbasis("test", l)

Test function. Usage: (H, dHx, dHy, dHz) = getSHbasis("test", 4)
"""
function getSHbasis(str::String, l::Int64)

	(nx,ny,nz) = (40,40,20)
	rx = range(-10,10,length=nx)
	ry = range(-10,10,length=ny)
	rz = range(-10,10,length=nz)
	(x,y,z) = ndgrid(rx,ry,rz)

	@time H = getSHbasis(x[:], y[:], z[:], l)
	@time (dHx, dHy, dHz) = getSHbasisGrad(x[:], y[:], z[:], l)

	nb = size(H,2)
	H = reshape(H, nx, ny, nz, nb)
	dHx = reshape(dHx, nx, ny, nz, nb)
	dHy = reshape(dHy, nx, ny, nz, nb)
	dHz = reshape(dHz, nx, ny, nz, nb)

	# jim(H[:,:,:,7], color=:jet)   # compare with >> evalspharm("test")

	return (H, dHx, dHy, dHz)
end

