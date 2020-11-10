using LinearAlgebra

"""
   shimmodel(r, A, s)

Implements f = HAs at spatial locations `r`, i.e., 
calculates the fieldmap (Hz) at `r` produced by shim settings `s`.

# Inputs
- `r`   [N 3]                  Spatial locations (cm) (distance from scanner iso-center)
- `A`   [9 9]                  Shim calibration matrix (includes the B0 term)
- `s`   length-9 vector        Shim settings (scanner-dependent hardware units)

# Output
- `f`   length-N vector        fieldmap (Hz) at locations `r`
"""
function shimmodel(
	r::Array{<:Real,2},
	A::Array{<:Real,2}, 
	s::Vector{<:Real}
	)

	N = size(r,1)

	(x,y,z) = (r[:,1], r[:,2], r[:,3])

	h = [ones(length(x)) x y z z.^2 x.*y z.*x x.^2-y.^2 z.*y];

	return [dot(h1, A, s) for h1 in eachrow(h)]
end
