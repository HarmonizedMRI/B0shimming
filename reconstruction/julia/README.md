# Julia B0 reconstruction

This folder reconstructs two interleaved Cartesian echoes from MATLAB-exported
k-space and estimates an unwrapped, regularized B0 field map.

## Input from MATLAB

Save the array currently named `d` as `kspace`, together with the echo times:

```matlab
kspace = single(d);
echo_times = single(TE);       % seconds
mask_threshold = single(0.4);  % optional

save('reconstruction_input.mat', ...
    'kspace', 'echo_times', 'mask_threshold', '-v7.3');
```

The expected k-space dimensions are:

```text
[nx, nEcho*ny, nz, ncoil]
```

Echoes are interleaved along the second dimension.

An optional 3D variable named `mask` can also be included. When it is absent,
a simple relative magnitude threshold is used.

## Run

### Interactive Julia session (recommended)

From this directory:

```bash
cd reconstruction/julia
julia
```

Activate the project and install the required packages (the first time only):

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Load the reconstruction module:

```julia
include("src/B0Reconstruction.jl")
using .B0Reconstruction
```

Estimate the field map:

```julia
estimate_fieldmap_file(
    "reconstruction_input.mat",
    "fieldmap.mat";
    fieldmap_nifti="fieldmap_hz.nii.gz",
    magnitude_nifti="magnitude_echo1.nii.gz",
)
```

### Command-line

Alternatively, the complete workflow can be run directly from the command line:

```bash
julia --project=. estimateFieldMap.jl \
    reconstruction_input.mat fieldmap.mat
```

The first run will automatically install the packages listed in `Project.toml`.

## Output

`fieldmap.mat` contains:

* `fieldmap_hz`: regularized field map from `MRIFieldmaps.jl`
* `unwrapped_fieldmap_hz`: ROMEO-unwrapped initialization
* `initial_fieldmap_hz`: wrapped initialization from `b0init`
* `magnitude_echo1`: first-echo magnitude image
* `mask`: processing mask
* `echo_times`: echo times in seconds

If the optional `fieldmap_nifti` and/or `magnitude_nifti` arguments are specified, the corresponding images are also written as NIfTI files (`.nii` or `.nii.gz`) using the acquisition geometry provided by the reconstruction.

## Notes

The current implementation:

* requires exactly two echoes;
* reconstructs each coil independently using a centered 3D inverse FFT;
* preserves the existing `shimtool.jl` complex-sum coil combination for field map estimation;
* uses a deliberately simple default mask that can later be replaced by the repository's application-specific masking workflow.

