# Code for estimating optimal B0 shim settings


## Using the Julia code: Quick start

1. Get this toolbox
```
$ git clone git@github.com:HarmonizedMRI/B0shimming.git
```

2. Change into the 'julia' subdirectory
```
cd julia
```

3. Start Julia (download from https://julialang.org/). Current version is 1.12.6.
```
julia
```
4. Press `]` to enter the Julia package manager and do:
```
(@v1.12) pkg> activate .
(julia) pkg> instantiate
```

5. Press `backspace` to get back to the Julia prompt.

6. If you are using the Siemens example, install the Python implementation of mapVBVD to load in twix files. In the terminal
```
python3 pip install pymapvbvd
```

7. Run the example:
```
julia> include("shimtool.jl")
```
  or
```
julia> include("brain_shim_siem.jl")
```
8. Answer 'y' through prompts

Each panel in the output image shows the field map (in Hz) before (left) and 
after (right) 2<sup>nd</sup> order shimming.
