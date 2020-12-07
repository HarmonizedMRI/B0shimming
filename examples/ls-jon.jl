#=
ls-adam.jl
solve LS problem using ADAM and Nesterov from Flux
=#

using LinearAlgebra: norm, opnorm
using Flux
using Random: seed!
using LaTeXStrings
using Plots; default(markerstrokecolor=:auto)

function ls_adam(y, A ;
	x0 = zeros(size(A,2)),
	opt = ADAM(), # default optimizer
	batchsize::Int = 1,
	niter::Int = 100,
	fun::Function = (x,iter) -> undef, # archiving function
)

	x = copy(x0)

	function shimmodel(A,x)
		return A*x
	end

	function	loss(y, A)
		return  1/2 * norm(y - shimmodel(A,x))^2 # LS cost function (no "x" arg!)
		# return 1/2 * norm(shimmodel(r,A,s) + f0)^2   
	end

	#loss = (y, A) -> 1/2 * norm(y - A*x)^2 # LS cost function (no "x" arg!)

	out = Array{Any}(undef, niter+1)
	out[1] = fun(x0, 0)
	θ = params(x) # magic here using objectid()
	data = [(y, A)] # passed as loss(data...) during optimization

	for iter = 1:niter
		Flux.train!(loss, θ, data, opt)
		out[iter+1] = fun(x, iter)
	end

	return x, out
end


if true # test
	seed!(0)
	M,N = 100, 80
	A = randn(M,N)
	L = opnorm(A)^2
	xtrue = randn(N)
	y = A * xtrue + randn(M)
	xh = A \ y # ideal LS solution
	x0 = zeros(N)

	cost = x -> 1/2 * norm(y - A*x)^2
	cmin = cost(xh)
	niter = 20
	fun = (x,iter) -> [cost(x) - cmin, norm(x-xh), time(), x]

	xa, oad = ls_adam(y, A ; x0=x0, niter=niter, fun=fun,
		#	opt = Descent(1/L),
			opt = ADAM(0.2), # hand-tuned learning rate
		)

	xn, ons = ls_adam(y, A ; x0=x0, niter=niter, fun=fun,
			opt = Nesterov(1/L), # requires Lipschitz constant
		)

	lf = x -> log10(max(x,1e-16)); grab = (o,i) -> hcat(o...)[i,:]
	lg = (o,i) -> lf.(grab(o,i))
	costk = out -> lg(out,1)
 	errk = out -> lg(out,2)
	timek = out -> grab(out,3) .- out[1][3]
 	allk = out -> (lg(out,1), lg(out,2), timek(out))
	cost_ad, err_ad, time_ad = allk(oad)
	cost_ns, err_ns, time_ns = allk(ons)

	plot(xlabel="iteration [k]", ylabel=L"\log_{10}(\Psi(x_k) - \Psi_*)")
	scatter!(0:niter, cost_ad, color=:red, label="cost adam")
	scatter!(0:niter, cost_ns, color=:blue, label="cost Nesterov")
	p1 = plot!(legend=:bottomleft)

	plot(xlabel="Time [s]", ylabel=L"\log_{10}(\|x_k - \hat{x}\|/\|\hat{x}\|)")
	scatter!(0:niter, err_ad, color=:red, label="NRMSD adam")
	scatter!(0:niter, err_ns, color=:blue, label="NRMSD Nesterov")
	p2 = plot!(legend=:bottomleft)

	display(plot(p1, p2))
end
