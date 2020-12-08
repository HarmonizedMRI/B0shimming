using NLopt

function loss(s::Vector, HA::Array, f0::Vector, shimlims::Tuple)
	
	(lin_max, hos_max, hos_sum_mxa) = shimlims

	shos = hos_max * atan.(s[5:9]) * 2/pi

	sa = vcat(s[1:4], shos)

   1/2 * norm(HA*sa + f0)^2
end


"""
s0: initial guess
HA: Basis H times calibration matrix A
f0: Acquired field map
shimlims[1]: max shim current on each linear shim coil
shimlims[2]: max shim current on each high-order shim coil
shimlims[3]: max combined high-order shim current
"""
function shimoptim(s0::Vector, HA::Array, f0::Vector, shimlims::Tuple)

	opt = Opt(:LD_MMA, length(s))
	opt.lower_bounds = vcat(-Inf, -shimlims[1]*ones(3,), -shimlims[2]*ones(5,))
	opt.upper_bounds = vcat(Inf,   shimlims[1]*ones(3,),  shimlims[2]*ones(5,))
	opt.xtol_rel = 1e-4

	opt.min_objective = loss

	(shat, out) = ls_adam(WHA, Wf0; s0=s0, opt=opt, niter=niter)

	shat

#	s = copy(s0)

#	ahos_max = 4000

#	function	loss(HA, f0)   # LS cost function (no "s" arg!)
#		shos = ahos_max * atan.(s[5:9]) * 2/pi
#		sa = vcat(s[1:4], shos)
#		1/2 * norm(HA*sa + f0)^2       
#	end

#	θ = params(s) # magic here using objectid()
#	data = [(HA,f0)] # passed as loss(data...) during optimization

#	out = Array{Any}(undef, niter+1)
#	out[1] = fun(s0, 0)

#	s[5:9] .= ahos_max * atan.(s[5:9]) * 2/pi

#	return s, out

end
