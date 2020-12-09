# Julia code for B0 shimming

# Example usage

1. Start Julia
2. Run demo script
   ```
   julia> include("shimdemo.jl");
	julia> shat                          # optimized shim currents
   ```

# Check SH basis

Compare the basis produced by getSHbasis.jl with basis obtained with two different Matlab implementations

```
>> evalspharm("test");   % plots bases obtained with evalspharm.m and poly_harmonic.m
julia> include("getSHbasis.jl"); H = getSHbasis("test");
julia> jim(H[:,:,:,8], color=:jet)     # (l,m) = (2,2) etc
```

Same shapes but normalized differently.
