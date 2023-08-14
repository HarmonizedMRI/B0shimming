# B0 shimming example from start to finish

[under construction]


## Overview

This folder contains a complete workflow for harmonized B0 field mapping. 
The workflow involves the following steps:

1. **Calibrate your scanner's B0 shimming channels**.
This is done by scanning a uniform phantom with a Pulseq sequence that we provide 
(see the [sequence/Pulseq](sequence/Pulseq) folder), 
and saving the acquired data to a file named *shimcal.mat*.
This just needs to be done once for each scanner.

2. **Acquire and reconstruct a B0 field map in the object you wish to shim over**.
This involves running a Pulseq scan and calculating the B0 field map,
and saving the field map to a file named *f0.mat*.

3. **Define the shim region**.
Here you define the region(s) (within the B0 field map `f0`) that you wish 
to shim over, and save the corresponding (binary) mask to a file named *shimvol.mat*.

4. **Calculate and apply the optimal shim current settings**.
This involves running *shim.jl* to obtain the new shim settings;
That script loads the *.mat* files created in the previous steps, 
and prints the reommended shim settings to the console.


The workflow just described involves creating the following experimental files:

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



## Create *shimcal.mat*

1. **Create the Pulseq sequence file (*b0.seq*).**
    This step involves executing the MATLAB script *writeB0.m*
    (in [../sequence/Pulseq/](../sequence/Pulseq/))
    to create the Pulseq file *b0.seq*.
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
    Run *b0.seq* multiple times using your vendor platform's Pulseq interpreter, 
    each with only one shim channel turned on.
    For each shim channel, acquire two B0 maps:
    first with amplitude `a/2`, then with `-a/2`. 
    `a` is chosen to avoid phase wraps in the individual B0 maps.
    This 'balanced' approach minimizes sensitivity to B0 drift.

    1. **GE users** may use the following settings:
        ```
        # GE
        a = 20       # x/y/z shims
        a = 1000     # 2nd order shims
        ```
        The corresponding matrix `S` is:
        ```
        # GE
        S = diag([<x> <y> <z> <z2> <xy> <zx> <x2y2> <zy>])
          = diag([20  20  20  1000 1000 1000 1000   1000])
        ```
        The script *shimcal_ge.pl* shows how the calibration data can be obtained on GE scanners,
        in an automated way and without having to manually set each shim channel amplitude.


    2. **Siemens users** may use the following settings:
        ```
        # Siemens
        a = 20       # x/y/z shims
        a = 200      # 2nd order shims
        ```
        The corresponding matrix `S` is:
        ```
        # Siemens
        S = diag([<x> <y> <z> <z2> <xy> <zx> <x2y2> <zy>])
          = diag([20  20  20  200  200  200  200    200 ])   
        ```

3. **Construct F and S, and write to file.**
    This involves reconstructing the (pairwise subtracted) B0 maps 
    and assembling the maps into the matrix `F`. 
    We also construct the shim amplitude matrix `S`, and define the object mask `mask_c`.

    1. **GE users** can perform these steps with the *makeshimcal_ge.m* script in this folder:
        ```
        >> makeshimcal_ge;    # Assumes P-files acquired with shimcal_ge.pl
        ```

    2. **Siemens users**...
       [TBD]
    

## Create *f0.mat*

1. Run *b0.seq* in the object we wish to shim over.
2. Reconstruct a B0 map and write to file:
    ```
    [TBD]
    ```
<!--
    >> getb0init;  % b0init, mask, magraw. Phase unwrapping is done in unwrap/main.jl
    << makef0;
-->


## Create *shimvol.mat*

```
[TBD]
```

<!--
    >> makeshimvol;  % uses FSL (bet) to do skull stripping
-->


## Calculate new shim settings

1. cd into the folder containing the Julia code.
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
    This script will load the various *.mat* files created above,
    and output the recommended shim settings.


