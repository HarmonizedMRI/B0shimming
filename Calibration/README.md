# B0 shimming Calibration



## Overview

This folder contains a complete workflow for harmonized B0 field mapping. 
The workflow involves the following steps:

1. **Calibrate your scanner's B0 shimming channels**.
This is done by scanning a uniform phantom with a Pulseq sequence that we provide 
(see the [sequence/Pulseq](sequence/Pulseq) folder), 
and saving the acquired data to a file named *shimcal.mat*.
This just needs to be done once for each scanner.


## 1. Calibrate and create *shimcal.mat*

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
    Run *b0.seq* 16 times using your vendor platform's Pulseq interpreter, 
    each with only one shim channel turned on.
    For each shim channel, acquire two B0 maps:
    first with amplitude `a/2`, then with `-a/2`. 
    `a` is chosen to avoid phase wraps in the individual B0 maps.
    This 'balanced' approach minimizes sensitivity to B0 drift.

    1. **GE users** may use the following settings:
        ```
        # GE
        a = 20       # x/y/z shims
        a = 100     # 2nd order shims
        ```
        The corresponding matrix `S` is:
        ```
        # GE
        S = diag([<x> <y> <z> <z2> <xy> <zx> <x2y2> <zy>])
          = diag([20  20  20  100  100  100  100   100])
        ```
        The script *shimcal_ge.pl* shows how the calibration data can be obtained on GE scanners,
        in an automated way and without having to manually set each shim channel amplitude. This script contains GE private commands, so it is available upon request.


    2. **Siemens users** may use the following settings:
        ```
        # Siemens
        a = 20       # x/y/z shims
        a = 100      # 2nd order shims
        ```
        The corresponding matrix `S` is:
        ```
        # Siemens
        S = diag([<A11> <B11> <A10> <A20> <A21> <B21> <A22> <B22>])
          = diag([20      20    20   100   100   100   100   100 ])   
        ```
       Currently, the calibration data has to be aquired manually on Siemens scanners.
       The 16 calibration measurements need to be aquired ordered in time, changing the shim channel settings from left to right.
       This means prior to the first measurement subtracting a = 20 [mikroT/m] from the first shim channel, then for the second measurement adding a = 20 [mikroT/m] to the first shim channel.
       Attention for the 7. measurement and all consequtive measurements a = 100! 
       These values are suggested as larger values can cause phase wraps in the calibration data.

       On Siemens only updating the shim settings leaves the center frequency in an unclear condition. 
       Run 'adjvalidate -fre -set "centerfreq"' after each shim update to avoid automatic center frequency manipulation by the system. The center frequency can be optained by running 'adjvalidate -fre -get' prior to the shim update. 

       : adjvalidate -fre -get
       -> "centerfreq"
       : adjvalidate -tra -get
       -> "transV"
       : adjvalidate -shim -get
       -> "current shim values"
       : adjvalidate -shim -set -mp ["current shim values"]+[-10 0 0 0 0 0 0 0]
       : adjvalidate -tra -set "transV"
       : adjvalidate -fre -set "centerfreq"
       

3. **Construct F and S, and write to file.**
    This involves reconstructing the (pairwise subtracted) B0 maps 
    and assembling the maps into the matrix `F`. 
    We also construct the shim amplitude matrix `S`, and define the object mask `mask_c`.

    1. **GE users** can perform these steps with the *makeshimcal_ge.m* script in this folder:
        ```
        >> makeshimcal_ge;    # Assumes P-files acquired with shimcal_ge.pl
        ```

    2. **Siemens users** can perform these steps with the *makeshimcal_siemens.m* script in this folder:
        ```
        >> makeshimcal_siemens;    # Assumes .dat-files acquired in order (see makeshimcalsiemens.m)
        ```
    3. Move shimcal.mat to a folder accessible by the shimming Julia script ("shimtoo.jl" or brain_shimsiem.jl").
    


