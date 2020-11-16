# An open toolbox for B0 shimming 

WIP

##  Goal

To provide an alternative to the scanner's built-in B0 shimming routine,
so that the linear and high-order B0 shims can be set according to well-defined 
(and potentially application-specific) critera.
For example, the user may want to:
1. Minimize root-mean-square (RMS) B0 inhomogeneity over a user-specified 3D subvolume.
1. Minimize the maximum through-voxel B0 gradient, to reduce signal loss in T2\*-weighted imaging.

To do this we define the system shim model
```
f(s) = H*A*s + f0         
f:  [N 1]              fieldmap (Hz), where N = number of voxels
f0: [N 1]              observed 'baseline' field map, e.g., after setting all shim currents to zero
H:  [N nb]             spherical harmonic basis (see getSHbasis.m) evaluated at the N voxel locations
A:  [nb nShim]          calibration matrix, where nb depends on SH degree
s:  [nShim 1]          change in shim current amplitudes from baseline (hardware units)
```

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
F: [N nShim]                          fieldmaps (Hz) obtained by turning on/off individual shim coils
S: [nShim nShim]                      applied shim currents (pairwise differences) used to obtain F
A = inv(H'*H)*H'*F*inv(S);            [nb nShim] 
```

Example:
```
>> % reconstruct unwrapped field maps F = [nx ny nz 8], obtained with shim settings S = [8 8]
>> nx = 64; ny = 64; nz = 64;                              % fieldmap matrix size
>> FOV = [20 20 20];                                       % fieldmap FOV (cm) 
>> [X,Y,Z] = getgrid(nx,ny,nz,FOV);                        % [nx ny nz], in same units as FOV
>> % define mask = [nx ny nz] (object support)
>> H = getSHbasis(X(mask),Y(mask),Z(mask),2);              % [N 9] where N = numel(X(mask))
>> % similarly, mask F and reshape to [N 8]
>> % form 8x8 matrix S containing (pairwise differences in) shim currents used to obtain F
>> A = getcalmatrix(F, H, S);
```


## How to perform 2nd order shimming

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

TODO
