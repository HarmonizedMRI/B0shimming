# getSHbasis.jl
using SphericalHarmonics: computeYlm
using CoordinateTransformations: SphericalFromCartesian
using MIRTjim: jim

# export numSH getSHbasis


"""
    c2sph(r, l, m)

Evaluate the `(l,m)` spherical harmonic at one spatial location `(x, y, z)`.
"""
function c2sph(x, y, z, l, m)
	a = SphericalFromCartesian()([x,y,z]) # todo: SVector?
	(rho, ϕ, θ) = (a.r, a.θ, π/2 - a.ϕ)   # (radius, azimuth, colatitude)
	Y = computeYlm(θ, ϕ; lmax=l)
	rho^l * Y[(l,m)]
end


"""
    n = numSH(L::Int)
Number of spherical harmonics for given `L`
i.e., `sum(2*(0:L) .+ 1)`, which is 9 for `L=2`.
Order of 0th-2nd order terms (for `L=2`):
* 1 cf (center frequency, Hz)
* 2 z
* 3 x
* 4 y
* 5 z2
* 6 zx
* 7 zy
* 8 x2y2
* 9 xy
"""
function numSH(L::Int)
	sum(2*(0:L) .+ 1)
end


"""
    getSHbasis!(h, x, y, z; L=2)

Store in `h` the spherical harmonic basis up to order `L` (default 2),
evaluated at spatial location `(x, y, z)`.
Vector `h` must have length `numSH(L)`.
"""
function getSHbasis!(h::AbstractVector{<:Real}, x::Real, y::Real, z::Real; L::Int = 2)

	ic = 1
	for l = 0:L, m = 0:l
		f = c2sph(x, y, z, l, m)
		h[ic] = real(f)
		ic += 1
		if m != 0
			h[ic] = imag(f)
			ic += 1
 		end
 	end

	h[1] = 1.0 # basis for center frequency offset

	return h
end


"""
    H = getSHbasis(x, y, z; L::Int=2)

Get spherical harmonic basis up to order `L` (default 2)
evaluated at spatial location vectors `x, y, z`
See `numSH` for ordering.
Output size is `(length(x), length(y), length(z), numSH(L))`
"""
function getSHbasis(
	x::AbstractVector{<:Real},
	y::AbstractVector{<:Real},
	z::AbstractVector{<:Real};
	L::Int=2
)

	T = promote_type(eltype.([x, y, z])..., Float32) # at least Float32
	H = zeros(T, length(x), length(y), length(z), numSH(L))
	for ix=1:length(x), iy=1:length(y), iz=1:length(z)
		getSHbasis!((@view H[ix,iy,iz,:]), x[ix], y[iy], z[iz]; L)
	end
	return H
end


# convenience method for scalar inputs
getSHbasis(x::Real, y::Real, z::Real; L::Int=2) = getSHbasis([x], [y], [z]; L)


"""
    H = getSHbasis("test"; L::Int = 2)

Test function
"""
function getSHbasis(str::String; L::Int = 2,
	nx::Int=22, ny::Int=20, nz::Int=18, fov::NTuple{3,Real}=(20,20,20),
)
	rx = LinRange(1,-1,nx)*fov[1]/2
	ry = LinRange(1,-1,ny)*fov[2]/2
	rz = LinRange(1,-1,nz)*fov[3]/2
	@time H = getSHbasis(rx, ry, rz; L)
	jim(H; ncol=numSH(L), color=:jet, gui=true, line3plot=false) # compare with >> evalspharm("test")
	return H
end

# getSHbasis("test"; L = 2);
