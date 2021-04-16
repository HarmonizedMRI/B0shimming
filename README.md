# An open toolbox for B0 shimming 

## Quick start

1. Start Julia (download from https://julialang.org/)
2. Press `]` to enter the Julia package manager and do:
```
(@v1.5) pkg> activate .
(@v1.5) pkg> instantiate
```
3. Press `backspace` to get back to the Julia prompt.
3. Run the demo program:
```
julia> cd("julia")
julia> include("shim.jl")
```


##  Goal

To provide an alternative to the scanner's built-in B0 shimming routine,
so that the linear and high-order B0 shims can be set according to well-defined 
(and potentially application-specific) critera.
For example, the user may want to minimize root-mean-square (RMS) B0 inhomogeneity 
over a user-specified 3D subvolume.

To do this we define the shim system model
```
f(s) = H*A*s + f0         
f:  [N 1]        fieldmap (Hz), where N = number of voxels
f0: [N 1]        observed 'baseline' field map, e.g., after setting all shim currents to zero
H:  [N nb]       spherical harmonic basis (see julia/getSHbasis.jl). nb = number of basis functions.
A:  [nb nb]      shim coil expansion coefficients for basis in H (see julia/getcalmatrix.jl)
s:  [nShim+1 1]  change in center frequency (cf) offset and shim current amplitudes from baseline (hardware units)
```
For 2nd order shim systems, nShim = 8 (3 linear and 5 2nd order).  
Each term in `H` is an `[N 1]` vector, evaluated at the same `N` spatial locations as `f`. 
The first column corresponds to the center frequency offset.

The goal here is to set the shim current vector `s` to make `f(s)` as homogeneous
as possible -- or more generally, to choose `s` according to some desired property of `f`
such as minimizing roughness or the maximum through-voxel gradient.

To do this we need to first **calibrate** the shim system to obtain `A`.


## Shim calibration (obtaining A)

We obtain `A` by turning the shims on/off one-by-one and acquiring a 3D fieldmap for each shim setting.
This can be done in a stationary phantom, and only needs to be done once for each scanner.
We then obtain `A` as follows:
```
F = HAS
F: [N nShim]                   fieldmaps (Hz) obtained by turning on/off individual shim coils
S: [nShim nShim]               applied shim currents (pairwise differences) used to obtain F
A = inv(H'*H)*H'*F*inv(S);     [nShim+1 nShim+1]  (see julia/getcalmatrix.jl). Include cf offset term.
```

See `julia/shim.jl` for a complete example, and additional information for how to construct F.



## How to perform 2nd order shimming

See also ./examples/demoWLS.m.

1. Acquire a baseline fieldmap `f0`. This does not need to use the same matrix size or FOV as the calibration scan.
2. Define a binary `mask` (control points), and reshape `f0` to `[N 1]` where `N = numel(X(mask))`
3. Construct the basis `H` on the control points, e.g.,
```
  >> nx = 64; ny = 64; nz = 20;                    % fieldmap matrix size
  >> FOV = [20 20 10];                             % fieldmap FOV (cm) 
  >> [X,Y,Z] = getgrid(nx,ny,nz,FOV);              % [nx ny nz], in same units as FOV
  >> H = getSHbasis(X(mask),Y(mask),Z(mask),2);    % [N 9] where N = numel(X(mask))
```
4. Load the calibration matrix `A`.
5. Apply the model `f = HAs + f0` and solve for `s` according to your desired loss function.
For example, to minimize the weighted RMS B0 field the Matlab solution is

```
  >> s = -(W*H*A)\(W*f0); 
```
where `W` is a diagonal spatial weighting matrix.


## Demo the code

```
>> cd examples;
>> shimdemoWLS;
```

