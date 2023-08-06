# B0 shimming example from start to finish

WIP

## Overview 

We first create the following experimental files:
```
<filename>      <variables>          <purpose>
shimcal.mat     F S FOV_c mask_c     Used to calculate the shim calibration matrix `A`
f0.mat          f0 FOV               Field map (Hz) and FOV (cm) we wish to shim over
shimvol.mat     mask                 Shim volume; logical/binary mask on grid defined by f0
```

The variables contained in these files are:
```
F        [N_c nShim]         Fieldmaps (Hz) obtained by turning on/off individual shim coils
                             N_c = number of voxels in calibration imaging volume (FOV_c)
                             nShim = number of shim channels (e.g., 3 (linear) or 8 (2nd order))
S        [nShim nShim]       Applied shim currents (pairwise differences) used to obtain F
FOV_c    [3]                 Field of view (cm) for the calibration data
mask_c   [nx_c ny_c nz_c]    Logical/binary mask indicating object support (N_c = nx_c*ny_c*nz_c)
f0       [nx ny nz]          Fieldmap we wish to shim over (Hz)
FOV      [3]                 FOV (cm) corresponding to f0
mask     [nx ny nz]          Logical/binary mask defining the desired shim region
```
Here, the subscript `_c` refers to the calibration data.

We then calculate the optimal shim settings by running the Julia
script ../julia/shim.jl.


## Create shimcal.mat

[TODO]

1. Run b0.seq with different shim settings
2. Assemble images into F, and the shim amplitudes used to collect F into S
3. Define `FOV_c` and `mask_c`.
4. Create shimcal.mat:
```
>> makeshimcal;
```

## Create f0.mat

[TODO]

1. Run b0.seq
2. Reconstruct the field map `f0`, define FOV, and write to f0.mat
```
>> getb0init;  % b0init, mask, magraw. Phase unwrapping is done in unwrap/main.jl
>> makef0;     % regularized B0 estimation done in Matlab or in b0reg/main.jl
```

## Create shimvol.mat 
```
>> makeshimvol;  % uses FSL (bet) to do skull stripping
```

## Calculate shims

1. cd into the folder containing the Julia code (main.jl)
1. Press `]` to enter the Julia package manager and do:
```
(@v1.7) pkg> activate .
(julia) pkg> instantiate
```
1. Press backspace to get back to the Julia prompt.
1. Run the script:
```
  julia> include("shim.jl")
```
This script will load the various .mat files created above,
and output the recommended shim settings.


