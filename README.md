# An open toolbox for vendor-neutral (harmonized) B0 shimming 

This repository aims to provide an alternative to the scanner's built-in B0 shimming routine,
so that the linear and high-order B0 shims can be set according to well-defined 
(and potentially application-specific) critera across different sites and vendor platforms.
We envision this tool as one component of a more harmonized cross-vendor MRI workflow 
in support of **reproducible MRI research**.

The key features of this toolbox are:

* **Vendor neutrality**: 
The entire workflow, from data acquisition to field map estimation, uses open-source and vendor-neutral tools
that are designed to ensure consistent and reproducible B0 shimming across sites and MRI scanner vendors.
Specifically, data acquisition is based on [Pulseq](https://pulseq.github.io/),
which guarantees that low-level details of the acquisition sequence
(sequence timing, RF/gradient pulse shapes, etc)
are *identical* across sites and vendor platforms.
In addition, all subsequent data processing is done using code provided in this repository,
in a completely vendor-independent way.

* **Robust field map estimation:** 
We estimate field maps in a robust way using the `b0map()` function from
[MRIFieldmaps.jl](https://github.com/MagneticResonanceImaging/MRIFieldmaps.jl).

* **Freedom to define the shim optimization critera**:
The framework allows for arbitrary (e.g., nonlinear) loss functions, 
and may be useful for exploring alternative shimming criteria (beyond least-squares) in the future. 
For example, the user may want to minimize root-mean-square (RMS) B0 inhomogeneity 
over a user-specified (and not necessarily contiguous) 3D subvolume.


## Overview and example usage

The B0 shimming workflow proposed here involves the following steps:

1. **Calibrating your scanner's B0 shimming channels**.
This is done by scanning a uniform phantom with a Pulseq sequence that we provide 
(see [sequence/Pulseq](sequence/Pulseq)). 
This just needs to be done once for each scanner.

2. **Acquiring and reconstructing a 3D B0 field map in the object you wish to shim over**.
This involves running a Pulseq scan and calculating the field map.

3. **Defining the shim region you wish to shim over**.
This can be done in a variety of ways, e.g., using the 
FSL plugin from https://shimming-toolbox.org/.

4. **Calculating and applying the optimal shim current settings**.
This involves running the Julia script *shim.jl* that we provide
(see the [julia](./julia) folder).
This script calculates optimal shim settings based on the provided
shim calibration data, 3D field map, and shim region.

To see precisely how these pieces fit together, see the [example](./example) folder.


## Code description

Understanding the following information is not strictly necessary 
in order to use this toolbox, but may be helpful for troubleshooting
or for those who wish to modify or contribute to this repository.

Here we use the following conventions (examples in parentheses):
* Variables are enclosed in back ticks in the Markdown code (`a`).
* File names are in italics (*shim.jl*).
* Scalars and vectors are in lower-case (`a`), while matrices are in upper-case (`F`).

The code in this repository is based on the following model:
```
f(s) = H*A*s + f0         
f:  [N]          fieldmap (Hz), where N = number of voxels
f0: [N]          observed 'baseline' field map, e.g., with your scanner's default shim setting.
H:  [N nb]       spherical harmonic basis (see julia/getSHbasis.jl). nb = # of basis functions.
A:  [nb nb]      shim channel expansion coefficients for basis in H (see julia/getcalmatrix.jl)
s:  [nShim+1]    change in center frequency (cf) and shim currents from baseline (hardware units)
```
For 2<sup>nd</sup> order shim systems, nShim = 8 (3 linear and 5 2<sup>nd</sup> order).  
Each column in `H` is an `N`-vector, evaluated at the same `N` spatial locations as `f`. 
The first column corresponds to the center frequency offset.
This toolbox provides support for spherical harmonic basis functions of arbitrary order
(see *julia/getSHbasis.jl*), but the code should work equally well with other bases.

The goal here is to set the shim current vector `s` to make `f(s)` as homogeneous
as possible -- or more generally, to choose `s` according to some desired property of `f`
such as minimizing roughness or the maximum through-voxel gradient.

To do this we need to first **calibrate** the shim system to obtain `A`,
which contains the basis expansion coefficients for a particular shim system.
We do this by turning the shims on/off one-by-one and acquiring a 3D fieldmap for each shim setting,
and assembling that data into a matrix `F`:
```
F: [N_c nShim]     fieldmaps (Hz) obtained by turning on/off individual shim channels 
S: [nShim nShim]   applied shim currents (pairwise differences) used to obtain F
```
`F` should be obtained in a stationary phantom, and only needs to be acquired once for each scanner.
We then obtain `A` by fitting each column in `F`
to the basis in `H`, using least-squares fitting (backslash in Julia); see *julia/getcalmatrix.jl*.

<!---
See `julia/example.jl` for a complete example, and additional information for how to construct F.
-->



## Cite

[Jon-Fredrik Nielsen, Berkin Bilgic, Jason P. Stockmann, Borjan Gagoski, 
Jr-Yuan George Chiou, Jr, Lipeng Ning, Yang Ji, Yogesh Rathi, 
Jeffrey A. Fessler, Douglas C. Noll, and Maxim Zaitsev.
An open toolbox for harmonized B0 shimming.
Proc. Intl. Soc. Mag. Reson. Med. 2021; 3772.](https://index.mirasmart.com/ISMRM2021/PDFfiles/3772.html)

