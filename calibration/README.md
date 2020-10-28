# Calibration scan and spherical harmonic basis fit 

##  Model
```
   W*f = W*H*A*s         
   W: diag_sp(mask(:))    Real-valued weighting matrix. mask = [nx ny nz]
   f: [nx*ny*nz 1]        fieldmap (Hz)
   H: [nx*ny*nz 9]        spherical harmonic basis. See getSHbasis.m.
   A: [9 9]               calibration matrix. See getcalmatrix.m.
   s: [9 1]               shim amplitudes, in the following order: 
                          (b0 x y z z2 xy zx x2-y2 zy). Hardware units, except b0 (Hz).
```

The goal here is to obtain A.

For an example, see 
./calMatrices/um3t,inside,01Nov2019/main.m


## Calibration scan

The model here is
```
   F = HAS
   F: [nx*ny*nz 8]    fieldmaps (Hz) obtained by turning on individual shim coils. See shimcal.pl.
   S: [9 9]           applied shim values (see shimcal.pl and getcalmatrix.m)
```

The goal of the calibration scan is to obtain F.

To do this, we turn on/off individual shims and acquire B0 map for each shim setting (use the TOPPE sequence in ../psd/Cartesian/).

Steps: 
1. Prescribe the TOPPE sequence in ../psd/Cartesian
    * Download
    * Auto prescan
1. Scan:
```
$ perl shimcal.pl
```


## Define ROI (mask)

Define ROI to shim over. You may want to use the GUI in ./GUI/, which creates the file 'ROI.h5'.

```
   mask: [nx ny nz]   Real-valued weighting mask
```

Example (Matlab):
```
>> roi = toppe.getroi('ROI.h5');
>> seq = getParams();
>> mask1 = roi2mask(roi,seq);
>> readoutFile = '../pulseSequence/Cartesian/readout.mod';
>> [~,imsos] = toppe.utils.recon3dft(pfile, 'readoutFile', readoutFile);
>> mask2 = imsos > 0.1*max(imsos(:));
>> mask = mask1.*mask2;
```


## Get calibration matrix 

1. Reconstruct UNWRAPPED fieldmaps 'F' ([nx ny nz 8])
1. Get calibration matrix:
```
>> addpath ..
>> seq = getparams();              % struct containing various acquisition/experimental parameters
>> A = getcalmatrix(F,seq,mask);
```

Calibration matrices are saved in ./calMatrices/


## How to perform 2nd order shimming using A

```
   W*f = W*H*A*s         
   W: diag_sp(mask(:))    Real-valued weighting matrix. mask = [nx ny nz]
```

1. Acquire fieldmap with ../psd/Cartesian/
1. Reconstruct UNWRAPPED fieldmap 'fmap' ([nx ny nz])
1. Calculate shim amplitudes and create the file 'newshims.atp':
```
>> s = getshims(fmap, A, seq, mask);
```
1. Copy 'newshims.atp' to scanner and apply shims.
```
$ atp newshims.atp 
```
1. Set scanner center frequency to the value of s(1).  **TODO:** add B0 offset to atp file
1. (optional) Acquire fieldmap again and compare with predicted field map (the figure produced by getshims).


## Test the scripts

```
>> addpath ..
>> seq = getparams();
>> A = getcalmatrix('test',seq)
>> s = getshims('test',A,seq)
```


