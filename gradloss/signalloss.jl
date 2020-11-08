using LinearAlgebra
using Plots

include("getexample.jl")

"""
   function signalloss(r, g, Δ, A, Δs, te)

Calculate relative/normalized signal loss due to
intra-voxel dephasing from B0 inhomogeneity.

# Inputs
- `r`   vector of length 3     Spatial position (cm) (distance from scanner iso-center)
- `g`   vector of length 3     Observed B0 field gradient (Hz/cm) at position `r` 
- `Δ`   vector of length 3     Voxel size (cm)
- `A`   [9 9]                  Shim calibration matrix (includes the B0 term)
- `Δs`  vector of length 9     Change in shim settings (from those used when observing `g`)
- `te`  scalar                 Time to echo (sec)
"""
function signalloss(
	r::Vector{<:Real}, 
	g::Vector{<:Real}, 
	Δ::Vector{<:Real}, 
	A::Array{<:Real,2}, 
	Δs::Vector{<:Real},
	te::Real)

	# Change in field gradient at position r due to applied shim changes Δs
	# The shim basis functions are h = [1 x y z z.^2 x.*y z.*x x.^2-y.^2 z.*y];
	(x,y,z) = r;
	dh = [0 1 0 0 0   y z  2*x 0;    # dh/dx
			0 0 1 0 0   x 0 -2*y z;    # dh/dy
			0 0 0 1 2*z 0 x  0   y]    # dh/dz
	dhAΔs = [dot(dh[1,:], A, Δs), dot(dh[2,:], A, Δs), dot(dh[3,:], A, Δs)]   # Hz/cm
	#print("Change in B0 gradient (Hz/cm) due to applied shims: $dhAΔs Hz/cm\n")

	# phase gradient at te (cycles/cm)
	gw = (dhAΔs + g)*te;        

	# return relative signal
	return sinc(dot(Δ, abs.(gw)))
end

"""
   function signalloss(r, g, Δ, A, Δs, te)

	Accept vector of voxel locations r (and their B0 gradients g) as input.
"""
function signalloss(
	r::Vector{Vector{Float64}},    # TODO: why doesn't <:Real work here?
	g::Vector{Vector{Float64}}, 
	Δ::Vector{<:Real}, 
	A::Array{<:Real,2}, 
	Δs::Vector{<:Real},
	te::Real)

	return map((r,g) -> signalloss(r,g,Δ,A,Δs,te), r, g)  

end

"""
	function signalloss(str::String)

Toy test function.
"""
function signalloss(str::String)

	N = 300;
	r = Vector{Vector{Float64}}(undef,N)
	g = Vector{Vector{Float64}}(undef,N)
	for i in 1:N
		r[i] = [2,2,2];                      # cm
		g[i] = 100*(i-N/2-0.5)/(N/2)*[1,1,1];         # Hz/cm
	end

	Δ  = 0.3*[1, 1, 1]     # cm

	A = collect(I(9)*1.)

	# h = [1   x   y   z   z.^2 x.*y z.*x x.^2-y.^2 z.*y]
	Δs =  [0., 0., 0,  0., 0,   0,   0,   0,        0   ]  # change in shim settings
	te = 25e-3     # echo time (sec)
	
	f = signalloss(r, g, Δ, A, Δs, te);

	df = Vector{Float64}(undef,N)
	for i in 1:N
		df[i] = g[i][1]
	end
	display(plot(df,f))

end

"""
   function cost(r, g, Δ, A, Δs, te, w)

Calculate weighted cost.
See signalloss.jl

# Inputs
- `r`   length-N vector of voxel positions (each entry is a length-3 vector)
- `g`   length-N vector of B0 field gradients (each entry is a length-3 vector)
- `Δ`, `A`, `Δs`, `te`: see signalloss.jl
- `w`   length-N vector of weights
"""
function cost(
	r::Vector{Vector{Float64}}, 
	g::Vector{Vector{Float64}}, 
	Δ::Vector{Float64}, 
	A::Array{Float64,2}, 
	Δs::Vector{Float64},
	te::Float64, 
	w::Vector{Float64})

	#f = map((r,g) -> signalloss(r,g,Δ,A,Δs,te), r, g)  
	f = signalloss(r,g,Δ,A,Δs,te)  

	return norm((1 .- f).*w)/sqrt(length(f))
end


"""
	cost(str::String)

Test function.
"""
function cost(str::String)

	(r,g,Δ,A,Δs,te) = getexample()

	N = length(r)    # number of voxels

	# spatial weights
	w = Vector{Float64}(undef,N)
	for i in 1:N
		w[i] = 1.; #i/N
	end

   @time f = cost(r, g, Δ, A, Δs, te, w)

	print("N: $N, f: $f \n")

end
