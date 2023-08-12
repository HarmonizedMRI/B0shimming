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
Here, the subscript `_c` denotes 'calibration'.

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



# OLD NOTES

### Scan, then compare pre/post fieldmap

Scan, save as `P,b0,post.7`. Then:
```
>> compareb0;
```

## Shim calibration files

### UM3TMR750

See ~/github/jfnielsen/scanLog/20220907_UM3TMR750_B0shim


## 

The shim calculation script is HarmonizedMRI/B0shimming/julia/shim.jl, which
requires access to the three .mat files listed above.

1. Create `shimcal.mat`  
Examples:  
20220907_UM3TMR750_B0shim/main2.m  
20220608_UM1.5TSola_shimcal/makeshimcal.m

2. Acquire B0 map and save to `f0.mat`  
Examples:  
20220907_UM3TMR750_B0shim/main3.m  
20221215_UM1.5TSola_shim/makef0.m

3. Create `shimvol.mat` containing shim volume mask  
   NB! Do this **after** creating `f0.mat` (since mask is used)
   1. Run localizer scan and convert to `Localizer.h5` for displaying in GUI:
   ```
   >> addpath github/toppeMRI/SlicePlanner/LocalizerScan
   >> pfile2hdf('P,localizer.7', 1, 'readout_localizer.mod');
   ```
   1. Start GUI and export ROI to ROI.h5, and create symbolic link to ROI.h5
   ```
   $ cd ~/github/toppeMRI/SlicePlanner/GUI
   $ ./start
   ...
   $ cd ~/shimtmpfiles/
   $ ln -s ~/github/toppeMRI/SlicePlanner/GUI/ROI.h5 .
   ```
   1. Convert ROI to image mask:
   ```
   >> roi = toppe.getroi('~/shimtmpfiles/ROI.h5', 1);
   >> load ~/shimtmpfiles/f0   % FOV, mask
   >> N = size(f0);
   >> masktmp = toppe.roi2mask(roi, N, FOV);
   >> mask = 1.0 * logical(mask & masktmp);
   >> save ~/shimtmpfiles/shimvol.mat mask
   ```

4. Calculate shim settings by running `shim.jl`:
   1. Change to B0shimming/julia directory
   ```
   $ cd ~/github/HarmonizedMRI/B0shimming/julia
   ```
   1. Start Julia and run the shim calculation:
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

### Example 1

1. In Matlab:
```
cd ~/github/jfnielsen/scanLog/20220907_UM3TMR750_B0shim/
main2;   % create ~/shimtmpfiles/shimcal.mat
main3;   % create ~/shimtmpfiles/f0.mat
main4;   % create ~/shimtmpfiles/Localizer.h5
```

1. Use GUI to define ROI to shim over

1. Create ~/shimtmpfiles/shimvol.mat 
```
main5;
```

1. Run `shim.jl`


### Example 2
See ~/github/jfnielsen/scanLog/20220912_UM3TMR750_B0shim



## In case it's needed: here are notes from ../20220909_UM3TMR750_B0shim/

In the following, the labels R/L, A/P, S/I refer to GE's direction for
Head First, Supine patient orientation, i.e., the 'patient coordinate system'.

1. Scan ACR phantom w/ localizer scan, created with:
```
>> main1;   % FOV = [24 24 24] cm; 2mm iso; interpolated to [240 240 240] matrix
```
TOPPE v4 settings:  
FOV centered on iso-center  
Freq. dir = R/L

1. Orient Java GUI viewports exactly as on GE MR750 console for
Head First, Supine patient orientation:  
```
Axial view:    Patient L = view right; Patient A = view top  
Sagittal view: Patient P = view right; Patient S = view top  
Coronal view:  Patient L = view right; Patient S = view top
```
DONE: this was done by editing pfile2hdf.m (SlicePlanner/Localizer)

1. Choose 'universal' coordinate system (UCS) for B0 shimming (and other things)
so the various pieces work nicely together, i.e., GUI, 
Matlab ROI to xyz conversion (SlicePlanner/Matlab/roi2xyz.m),
Julia shim tool, etc. 
Let's go with a right-handed system with mm units:      
```
Patient R = +x; A = +y; S = +z; with (0,0,0) at iso-center  
unit: mm
```

1. Make ROI export and console dump UCS compliant -- DONE

1. Make scroll wheel in GUI move to higher (more Superior) slices when scrolling away/up -- DONE

1. Orient image volume (recon'd w/ toppe.utils.recon3dft.m) with UCS, so it appears correctly in 'im'  -- DONE (see main3.m)
Accordingly, added 'alignWithUCS' option to toppe.utils.recon3dft()
