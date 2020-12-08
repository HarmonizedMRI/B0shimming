using LinearAlgebra: norm, opnorm
using Flux
using MIRT: embed!
using Random: seed!
using LaTeXStrings
using Plots; default(markerstrokecolor=:auto)

function ls_adam(
	HA::Array{<:Real,2},
	f0::Vector{<:Real}; 
	s0 = zeros(size(HA,2)),
	opt = ADAM(), # default optimizer
	batchsize::Int = 1,
	niter::Int = 100,
	fun::Function = (s,iter) -> undef, # archiving function
)

	s = copy(s0)

	ahos_max = 4000

	function	loss(HA, f0)   # LS cost function (no "s" arg!)
		# norm(f0 + HA*s,6)^6
		shos = ahos_max * atan.(s[5:9]) * 2/pi
		sa = vcat(s[1:4], shos)
		1/2 * norm(HA*sa + f0)^2       
	end

	θ = params(s) # magic here using objectid()
	data = [(HA,f0)] # passed as loss(data...) during optimization

	out = Array{Any}(undef, niter+1)
	out[1] = fun(s0, 0)

	for iter = 1:niter
		Flux.train!(loss, θ, data, opt)
		out[iter+1] = fun(s, iter)
	end

	s[5:9] .= ahos_max * atan.(s[5:9]) * 2/pi

	return s, out
end

#	s0 = zeros(9)
#	cost = s -> 1/2 * norm(HA*s .+ f0)^2 
#	fun = (s,iter) -> cost(s) #, time()]
#	opt = ADAM(0.2)
#	@time (shat, out) = ls_adam(HA,f0; s0=s0, niter=1000, fun=fun, opt=opt)
#	matwrite("results.mat", Dict([("shat", shat), ("H", H), ("A", A), ("mask", mask), ("f0", f0)]))
