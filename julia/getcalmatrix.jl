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
	inv(H'*H)*H'*F*inv(S)
end
