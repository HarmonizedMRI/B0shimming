n = 500
p = 1000
x = randn(n,p)'
y = sum(x[1:5,:], dims=1) .+ randn(n)'*0.1

using AutoGrad
v = [0.0001*randn(1,p) 0.0]

loss2(v) = sumabs2(y - (v[1]*x .+ v[2])) / size(y,2)
lossgradient = grad(loss2)
function train(v, x, y, lambda; lr=.1)
        g = lossgradient(v)
        v[1] -= lr * g[1]
        v[2] -= lr * g[2]
        v[1][-(lr * lambda) .< v[1] .< (lr * lambda)] = 0
    return v
end

niter = 50
lambda = 1.0
@time for i=1:niter; v = train(v, x, y, lambda); println(loss2(v)); end

