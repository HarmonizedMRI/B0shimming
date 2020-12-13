# Julia code for B0 shimming

## Usage

1. Save calibration data and baseline fieldmap to .jld2 files, see `mat2jld2.jl` for an example.
1. Edit the top section in `shim.jl` accordingly
1. Calculate optimized shim settings
  julia> include("shim.jl")



## Check SH basis

Compare the basis produced by getSHbasis.jl with basis obtained with two different Matlab implementations

```
>> evalspharm("test");   % plots bases obtained with evalspharm.m and poly_harmonic.m
julia> include("getSHbasis.jl"); H = getSHbasis("test");
julia> jim(H[:,:,:,8], color=:jet)     # (l,m) = (2,2) etc
```

Same shapes but normalized differently.
