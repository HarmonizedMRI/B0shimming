using LinearAlgebra

"""
   function signalloss(r, g, Δ, A, Δs)

Calculate (relative/normalized) signal loss due to
intra-voxel dephasing from B0 inhomogeneity.

# Inputs
- `r`       [3 1]      Spatial position (cm) (distance from scanner iso-center)
- `g`       [3 1]      Observed B0 field gradient (Hz/cm) at position `r` 
- `Δ`       [3 1]      Voxel size (cm)
- `A`       [9 9]      Shim calibration matrix (includes the B0 term)
- `Δs`      [9 1]      Change in shim settings (from those used to acquire `g`)
"""
function signalloss(r::Vector{<:Real}, g::Vector{<:Real}, Δ::Vector{<:Real}, A::Array{<:Real,2}, Δs::Vector{<:Real})

	# Spherical harmonic basis evaluated at r
	(x,y,z) = r;
	h = [1 x y z z.^2 x.*y z.*x x.^2-y.^2 z.*y];

	# Change in field gradient due to applied shim changes Δs
	dh = [0 1 0 0 0 y z 2*x 0;    # dh/dx
			0 0 1 0 0 x 0 -2*y z;   # dh/dy
			0 0 0 1 2*z 0 x 0 y];   # dy/dz
	AΔs = A*Δs;
	dhAΔs = [dot(dh[1,:], AΔs), dot(dh[2,:], AΔs), dot(dh[3,:], AΔs)]
	# print("Additional B0 gradient due to applied shim changes: $dhAΔs\n");

	# return relative signal loss
	# print("Δ: $Δ, g: $g \n");
	return sinc(dot(Δ, dhAΔs + g));

end

# test function
function signalloss(str::String)

	r = [2,3,4];            # cm
	g = 0*[1,2,3];          # Hz/cm
	Δ  = [0.3, 0.3, 0.3];   # cm

	A = [1 0.1 0 -0.1 2.3 0 0 0 0;
		  0 1 0.001 0.002 0 0 0 0 0;
		  0 0 1 0 0 0 0 0 0;
		  0 0 0 1 0 0 0 0 0;
		  0 0 0 0 1 0 0 0 0;
		  0 0 0 0 0 1 0 0 0;
		  0 0 0 0 0 0 1 0 0;
		  0 0 0 0 0 0 0 1 0;
		  0 0 0 0 0 0 0 0 1];

	Δs = [0,-2,0,0,0,0,-0.2,0,0.2];  # change in shim settings
	
	f = signalloss(r, g, Δ, A, Δs);
	print("Relative signal: $f");
end
