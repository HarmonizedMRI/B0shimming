using MIRT: ndgrid, jim
using JLD2
using Plots

include("getSHbasis.jl")
include("fieldgrad.jl")

"""
	getSHbasisGrad(x, y, z, L)

Get gradient of spherical harmonic basis up to order L evaluated at spatial locations x, y, z.
Uses one-sided finite differences (faster than ForwardDiff.gradient)
"""
function getSHbasisGrad(
	x::Vector{<:Real}, 
	y::Vector{<:Real}, 
	z::Vector{<:Real}, 
	L::Int64
	)

	d = 1e-3    # distance step

	H = getSHbasis(x, y, z, L);
	dHx = getSHbasis(x .+ d, y,      z,      L) - H;
	dHy = getSHbasis(x,      y .+ d, z,      L) - H;
	dHz = getSHbasis(x,      y,      z .+ d, L) - H;
	
	gHx = dHx/d
	gHy = dHy/d
	gHz = dHz/d

	return (gHx, gHy, gHz)
end

"""
	getSHbasisGradAuto(x, y, z, L)
"""
function getSHbasisGradAuto(
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
	getSHbasisGrad("test", l::Int64)

Test function.   
Compares `gHx*A*s` with gradient of `H*A*s` -- will differ slightly due to finite grid size. 

"""
function getSHbasisGrad(str::String, l::Int64)

	(nx,ny,nz) = (30,30,30)
	fov = [10, 10, 10]             # cm

	(x,y,z) = ndgrid(
		range(-fov[1]/2, fov[1]/2, length=nx), 
		range(-fov[2]/2, fov[2]/2, length=ny), 
		range(-fov[3]/2, fov[3]/2, length=nz) 
	)

	@time (gHx, gHy, gHz) = getSHbasisGrad(x[:], y[:], z[:], l)
	#@time (gHxa, gHya, gHza) = getSHbasisGradAuto(x[:], y[:], z[:], l)

	@time H = getSHbasis(x[:], y[:], z[:], l)

	@load "A_l2.jld2"    # see shim.jl

	# shim amplitudes
	s = [30, 4, 2, 10, -1732, 660, -83, 1646, -230]
	#s = [30; rand(-10:10,3); rand(-2000:2000, 5)]
	@show s

	# Calculate gradient in two ways and compare

	# Calculate fieldmap produced by shims 's' and take finite diff of result
	f = H*A*s
	f = reshape(f, nx, ny, nz)
	@time (gx, gy, gz) = fieldgrad(fov, f)

	# Calculate fieldmap using the gradients of the basis H
	nb = size(gHx,2)
	gx2 = gHx*A*s
	gy2 = gHy*A*s
	gz2 = gHz*A*s
	gx2 = reshape(gx2, nx, ny, nz)
	gy2 = reshape(gy2, nx, ny, nz)
	gz2 = reshape(gz2, nx, ny, nz)

	clim = (-5,5)
	p1 = jim(cat(gx-gx2;dims=1); clim=clim, color=:jet)
	p2 = jim(cat(gy-gy2;dims=1); clim=clim, color=:jet)
	p3 = jim(cat(gz-gz2;dims=1); clim=clim, color=:jet)

	p = plot(p1, p2, p3, layout=(3,1))
	pyplot()
	display(p)

	return (gx,gy,gz,gx2,gy2,gz2)
end

