
function shimoptim(s0::Vector, WHA::Array, Wf0::Vector) # , max_lin::Real, max_hos::Real, max_sum_hos::Real)

	opt = ADAM(0.4)
	niter = 300
	ahos_max = 4000

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
