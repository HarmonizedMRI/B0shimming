x = Param([1,2,3])
y = @diff sum(abs2,x)
grad(y,x)



n = 500
p = 1000
x = randn(n,p)'
y = sum(x[1:5,:], dims=1) .+ randn(n)'*0.1

using AutoGrad
w = [0.0001*randn(1,p) 0.0]

loss(w) = sumabs2(y - (w[1]*x .+ w[2])) / size(y,2)
lossgradient = grad(loss)
function train(w, x, y, lambda; lr=.1)
        g = lossgradient(w)
        w[1] -= lr * g[1]
        w[2] -= lr * g[2]
        w[1][-(lr * lambda) .< w[1] .< (lr * lambda)] = 0
    return w
end

niter = 50
lambda = 1.0
@time for i=1:niter; w = train(w, x, y, lambda); println(loss(w)); end

