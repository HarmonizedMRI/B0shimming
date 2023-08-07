# B0 shimming example from start to finish

[under construction. Jon 6-Aug-2023]

TODO:
* Streamline workflow by calling Julia from MATLAB? Looks like not so straightforward.
* How to streamline Julia unwrap step (awkward to go back to MATLAB before running shim.jl)
* Creating f0.mat: Update scripts, and make them vendor-agnostic
   * getb0init.m: clean up, and separate out Julia unwrap step
   * makef0.m  (and separate out regularized estimation part)
* update makeshimvol.m
* Comment on how to define shim region.


## Overview 

1. We first create the following experimental files:
    ```
    <filename>      <variables>          <purpose>
    shimcal.mat     F S FOV_c mask_c     Used to calculate the shim calibration matrix `A`
    f0.mat          f0 FOV               Field map (Hz) and FOV (cm) we wish to shim over
    shimvol.mat     mask                 Shim volume; logical/binary mask on grid defined by f0
    ```

    The variables contained in these files are:
    ```
    F   [nx_x ny_c nz_c nShim]  Fieldmaps (Hz) obtained by turning on/off individual shim coils
                                N_c = number of voxels in calibration imaging volume (FOV_c)
                                nShim = number of shim channels (e.g., 3 (linear) or 8 (2nd order))
    S        [nShim nShim]      Applied shim currents (pairwise differences) used to obtain F
    FOV_c    [3]                Field of view (cm) for the calibration data
    mask_c   [nx_c ny_c nz_c]   Logical/binary mask indicating object support (N_c = nx_c*ny_c*nz_c)
    f0       [nx ny nz]         Fieldmap we wish to shim over (Hz)
    FOV      [3]                FOV (cm) corresponding to f0
    mask     [nx ny nz]         Logical/binary mask defining the desired shim region
    ```
    Here, the subscript `_c` refers to the calibration data.

2. We then calculate the optimal shim settings by running the Julia
script ../julia/shim.jl.


## Create shimcal.mat

1. **Create the Pulseq sequence file.**
    This step involves executing the MATLAB script writeB0.m 
    (in [../sequence/Pulseq/](../sequence/Pulseq/))
    to create the Pulseq file `b0.seq`.
    For this you will need the Pulseq toolbox:
    ```
    $ git clone git@github.com:pulseq/pulseq.git
    ```
    Then in MATLAB, do:
    ```
    >> addpath pulseq/matlab
    >> writeB0;
    ```
2. **Acquire the data.**
    Run b0.seq multiple times using your vendor platform's Pulseq interpreter, 
    each with only one shim channel turned on.
    For example, for GE scanners, the following settings may be used:
    ```
    # GE
    S = diag([<x> <y> <z> <z2> <xy> <zx> <x2y2> <zy>])
      = diag([20  20  20  1000 1000 1000 1000   1000])
    ```
    where each entry denotes the shim current amplitude for a particular channel.
    These settings are chosen so as to avoid phase wraps in the individual B0 maps.
    For Siemens, the following settings may be used:
    ```
    # Siemens
    S = diag([<x> <y> <z> <z2> <xy> <zx> <x2y2> <zy>])
      = diag([20  20  20  200  200  200  200    200 ])   
    ```
    In addition, we acquire one **reference image volume** with all shim channels off.

    [NB! When testing we actually acquired two B0 maps for each shim channel:
    first with amplitude 'a/2', then with '-a/2'. 
    The shim amplitudes listed above are the **difference** between these settings.
    Let's refer to this as a 'balanced' acquisition. Not sure if we need to do this?]
    
3. **Construct F and S, and write to file.**
    This involves reconstructing the B0 maps 
    (after subtracting the phase in the reference acquisition)
    and assembling the maps into the matrix `F`. 
    We also construct the shim amplitude matrix `S`, and define the object mask `mask_c`.
    In this example, we perform these steps with the 'makeshimcal.m' script in this folder:
    ```
    >> makeshimcal;
    ```
    

## Create f0.mat

1. Run b0.seq in the object we wish to shim over.
2. Reconstruct a B0 map and write to file:
```
>> getb0init;  % b0init, mask, magraw. Phase unwrapping is done in unwrap/main.jl
>> makef0;     % regularized B0 estimation done in Matlab or in b0reg/main.jl
```


## Create shimvol.mat 
```
>> makeshimvol;  % uses FSL (bet) to do skull stripping
```


## Calculate new shim settings

1. cd into the folder containing the Julia code (main.jl)
1. Make sure the files you created above are accessible from this folder, 
   e.g., create copies or symbolic links:
   ```
   $ ln -s ~/mydata/shimcal.mat .
   $ ln -s ~/mydata/f0.mat .
   $ ln -s ~/mydata/shimvol.mat .
   ```
1. Start a Julia session 
    ```
    $ julia
    ```
1. Press `]` to enter the Julia package manager, then do:
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


