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

	return norm((1 .- f).*w)/length(f)
end


"""
	cost(str::String)

Test function.
"""
function cost(str::String)

	(Δ,A,Δs,te)  = getexample()

	# spatial weights
	w = Vector{Float64}(undef,N)
	for i in 1:N
		w[i] = 1.; #i/N
	end

   @time f = cost(r, g, Δ, A, Δs, te, w)

	print("N: $N, f: $f \n")

end

