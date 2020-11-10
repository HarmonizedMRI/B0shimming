using LinearAlgebra

"""
   shimmodel(r, A, s)

Implements f = HAs at spatial location `r`, i.e., 
calculates the fieldmap (Hz) at `r` produced by shim settings `s`.

This function is overloaded and also accepts a vector of spatial location vectors.

# Inputs
- `r`   vector of length 3     Spatial position (cm) (distance from scanner iso-center)
- `A`   [9 9]                  Shim calibration matrix (includes the B0 term)
- `s`   vector of length 9     Shim settings (scanner-dependent hardware units)

# Output
- `f`   scalar                 fieldmap (Hz) at spatial location `r`
"""
function shimmodel(
	r::Vector{<:Real}, 
	A::Array{<:Real,2}, 
	s::Vector{<:Real})

	(x,y,z) = r;
	h = [1 x y z z.^2 x.*y z.*x x.^2-y.^2 z.*y];

	return dot(h, A, s)   
end


"""
   shimmodel(r, A, s)

Acceps a vector of location vectors.
"""
function shimmodel(
	r::Vector{Vector{Float64}},    # TODO: why doesn't <:Real work here?
	A::Array{<:Real,2}, 
	s::Vector{<:Real})

	return map((r) -> shimmodel(r,A,s), r)  
end
