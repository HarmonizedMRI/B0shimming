using NLopt

"""
HA: Spherical harmonic basis H times calibration matrix A
f0: Acquired field map
shimlims[1]: max shim current on each linear shim coil
shimlims[2]: max shim current on each high-order shim coil
shimlims[3]: max combined high-order shim current
s0: initial guess
"""
function shimoptim(HA::Array, f0::Vector, shimlims::Tuple; s0::Vector=zeros(size(HA,2),))

	(lin_max, hos_max, hos_sum_max) = shimlims

	function loss(s::Vector, grad::Vector, HA, f0)
   	1/2 * norm(HA*s + f0)^2
	end

	opt = Opt(:LN_COBYLA, length(s0))
	opt.lower_bounds = vcat(-Inf, -lin_max*ones(3,), -hos_max*ones(5,))
	opt.upper_bounds = vcat( Inf,  lin_max*ones(3,),  hos_max*ones(5,))
	opt.xtol_rel = 1e-4
	opt.min_objective = (x, grad) -> loss(x, grad, HA, f0)

	(minf,minx,ret) = optimize(opt, s0)

	numevals = opt.numevals # the number of function evaluations
	println("got $minf at $minx after $numevals iterations (returned $ret)")

	s = minx
end
