# Run calibration scan for GE

The model here is
```
   F = CAS
   F: [nx*ny*nz 8]    fieldmaps (Hz) obtained by turning on individual shim coils. See shimcal.pl.
   S: [9 9]           applied shim values (see shimcal.pl and getCalMatrix.m)
```

The goal of the calibration scan is to obtain F.

To do this, we turn on/off individual shims and acquire B0 map for each shim setting (use the sequence in /..) 

Steps:
1. Prescribe the TOPPE sequence 
    * Download
    * Auto prescan
1. Scan:
```
$ perl shimcal.pl
```

(Contact jfnielse@umich.edu for a copy of shimcal.pl)

