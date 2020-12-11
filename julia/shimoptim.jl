using NLopt

"""
	Optimize shims using constrained nonlinear optimization and arbitrary loss function

function shimoptim(HA::Array, f0::Vector, shimlims::Tuple; s0::Vector, loss::Function)

# Inputs
- `HA`: Spherical harmonic basis H times calibration matrix A  
- `f0`: Acquired field map
- `shimlims`
  - shimlims[1]: max shim current on each linear shim coil
  - shimlims[2]: max shim current on each high-order shim coil
  - shimlims[3]: max combined high-order shim current
- `s0`: initial guess
- `loss`: loss function

# Output
- `s`: optimized shim currents (including DC offset)
"""
function shimoptim(HA::Array, f0::Vector, shimlims::Tuple; 
	s0::Vector = zeros(size(HA,2),),
	loss::Function = (s, HA, f0) -> 1/2 * norm(HA*s + f0)^2,
	ftol_rel = 1e-4
	)

	(lin_max, hos_max, hos_sum_max) = shimlims

	opt = Opt(:LN_COBYLA, length(s0))
	opt.ftol_rel = ftol_rel
	opt.min_objective = (s, grad) -> loss(s, HA, f0)

	opt.lower_bounds = [-Inf; -lin_max*ones(3,); -hos_max*ones(5,)]
	opt.upper_bounds = [ Inf;  lin_max*ones(3,);  hos_max*ones(5,)]
	opt.inequality_constraint = (s, grad) -> sum(abs.(s[5:9])) - hos_sum_max

	(minf,mins,ret) = optimize(opt, s0)

	#numevals = opt.numevals # the number of function evaluations
	#println("Loss=$minf after $numevals iterations (returned $ret)")

	mins
end
