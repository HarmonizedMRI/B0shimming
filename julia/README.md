# Julia code for B0 shimming

## Usage

1. Save calibration data and baseline fieldmap to .jld2 files, see `mat2jld2.jl` for an example.
1. Edit the top section in `shim.jl` accordingly
1. Calculate optimized shim settings
```
  julia> include("shim.jl")
```

Summary of steps:

```
julia> f0 = f0[mask]                                           # baseline B0 map before shimming
julia> l = 6                                                   # spherical harmonic order
julia> H = getSHbasis(x[mask], y[mask], z[mask], l)            # spherical harmonics evaluated at (x,y,z)
julia> @load "A.jld2"                                          # calibration matrix (see getcalmatrix.jl)
julia> max_lin = 100                                           # max linear shim current
julia> max_hos = 4000                                          # max high order shim current
julia> maxSumHOS = 12000                                       # max total HOS current
julia> lims = (max_lin*ones(3,), max_hos*ones(5,), maxSumHOS)
julia> lossfun = (s, HA, f0) -> 1/2*norm(HA*s + f0)^2          # loss function
julia> s = shimoptim(H*A, f0, lims; loss=lossfun)              # returns optimized shims 
```



## Check SH basis

Compare the basis produced by getSHbasis.jl with basis obtained with two different Matlab implementations

```
>> evalspharm("test");   % plots bases obtained with evalspharm.m and poly_harmonic.m
julia> include("getSHbasis.jl"); H = getSHbasis("test");
julia> jim(H[:,:,:,8], color=:jet)     # (l,m) = (2,2) etc
```

Same shapes but normalized differently.
