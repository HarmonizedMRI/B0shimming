using NLopt

"""
	Optimize shims using constrained nonlinear optimization and arbitrary loss function

function shimoptim(HA::Array, f0::Vector, shimlims::Tuple; s0::Vector, loss::Function)

# Inputs
- `HA`: Spherical harmonic basis H times calibration matrix A  
- `f0`: Acquired field map
- `shimlims`
  - shimlims[1]: (vector) max shim current on each linear shim coil
  - shimlims[2]: (vector) max shim current on high-order shim coils
  - shimlims[3]: (scalar) max combined high-order shim current
- `s0`: initial guess
- `loss`: loss function

# Output
- `s`: optimized shim currents (including center frequency offset)
"""
function shimoptim(HA::Array, f0::Vector, shimlims::Tuple; 
	s0::Vector = zeros(size(HA,2),),
	loss::Function = (s, HA, f0) -> 1/2 * norm(HA*s + f0)^2,
	ftol_rel = 1e-3,
	cflim = 2000    # empirical observation: needs to be finite for numerical stability
	)

	(lin_max, hos_max, hos_sum_max) = shimlims

	opt = Opt(:LN_COBYLA, length(s0))
	opt.ftol_rel = ftol_rel
	opt.min_objective = (s, grad) -> loss(s, HA, f0)

	# shim current limits
	if length(s0) > 4
		opt.lower_bounds = [-cflim; -lin_max; -hos_max]
		opt.upper_bounds = [ cflim;  lin_max;  hos_max]
	else
		opt.lower_bounds = [-cflim; -lin_max]
		opt.upper_bounds = [ cflim;  lin_max]
	end

	# limit on total HO shim current
	if ~isempty(hos_sum_max) & length(s0) > 4
		opt.inequality_constraint = (s, grad) -> sum(abs.(s[5:end])) - hos_sum_max
	end

	(minf,mins,ret) = optimize(opt, s0)

	#numevals = opt.numevals # the number of function evaluations
	#println("Loss=$minf after $numevals iterations (returned $ret)")

	mins
end
