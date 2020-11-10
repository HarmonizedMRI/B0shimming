# shimdemoFluxWLS.jl
#
# Same goal as shimdemoWLS.m, but solving for the shims using Flux in Julia.
# This general framework can then be used with other cost functions.

using LinearAlgebra: norm
using MIRT: embed!
using Flux

include("shimmodel.jl");

# get shim model and example data 
include("loadexampledata.jl");   # load A, f0, X, Y, Z, mask

# reformat/reshape data as needed 

opt = ADAM(0.2), # hand-tuned learning rate

function ls_adam(r, A, f0;
	s0 = zeros(size(A,2)),
	opt = ADAM(), # default optimizer
	batchsize::Int = 1,
	niter::Int = 100,
	fun::Function = (x,iter) -> undef, # archiving function
)

	s = copy(s0)
	loss = (r, A, f0) -> 1/2 * norm(shimmodel(r,A,s) + f0)^2   # LS cost function (no "s" arg!)
	out = Array{Any}(undef, niter+1)
	out[1] = fun(s0, 0)
	θ = params(s) # magic here using objectid()
	data = [(r,A,f0)]   # passed as loss(data...) during optimization

	for iter = 1:niter
		Flux.train!(loss, θ, data, opt)
		out[iter+1] = fun(s, iter)
	end

	return s, out
end



