# Code for estimating optimal B0 shim settings

[more here]

## Using the Julia code: Quick start

1. Get this toolbox
```
$ git clone git@github.com:HarmonizedMRI/B0shimming.git
```

2. Change into the 'julia' subdirectory
```
julia> cd("julia")
```

3. Start Julia (download from https://julialang.org/). Current version is 1.6.0.

4. Press `]` to enter the Julia package manager and do:
```
(@v1.6) pkg> activate .
(julia) pkg> instantiate
```

5. Press `backspace` to get back to the Julia prompt.

6. Run the example:
```
julia> include("example.jl")
```

Each panel in the output image shows the field map (in Hz) before (left) and 
after (right) 2<sup>nd</sup> order shimming of a cylindrical jar phantom:  
![output of examples.jl](resources/out30.png "Example")
