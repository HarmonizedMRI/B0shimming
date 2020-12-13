using LinearAlgebra

function getcalmatrix(
	F::Array{<:Real,2},
	H::Array{<:Real,2},
	S::Array{<:Real,2}
	)

	# Insert DC terms
	N = size(F,1)
	F = [ones(N) F] 
	Sin = S;
	S = zeros(9,9);
	S[1,1] = 1;
	S[2:end,2:end] = Sin

	# return calibration matrix
	lam = 1e2;
	A = inv(H'*H + lam*I(size(H,2)))*H'*F*inv(S)

	# center frequency offset basis is 1 by definition
	A[1,1] = 1.0;
	A[2:end,1] .= 0;

	return A
end
