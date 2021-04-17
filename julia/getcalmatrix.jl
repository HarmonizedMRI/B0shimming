using LinearAlgebra

"""
	getcalmatrix(F::Array, H::Array, s::Vector)

TODO: complete help text

"""
function getcalmatrix(
	F::Array{<:Real,2},   # [nvoxels nshims] Acquired field maps
	H::Array{<:Real,2},   # [nvoxels nbasis] Spatial basis functions
	s::Vector             # applied shim currents (e.g., pairwise differences) used to acquire F
	)

	nb = size(H,2)      # number of basis functions. Includes center frequency.
	ns = length(s)      # number of shim terms (excluding center frequency)

	A = zeros(nb, ns+1)

	# center frequency (identity by definition)
	A[1,1] = 1.0

	# fit each acquired fieldmap to basis
	for ii = 1:ns
		A[:,ii+1] = s[ii] * H \ F[:,ii]
	end

	# A = inv(H'*H)*H'*F*inv(S);   [nb nb] Includes cf offset term.	

	return A
end
