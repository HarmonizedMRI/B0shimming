# B0 field map reconstruction

This folder contains the reconstruction pipeline for the vendor-neutral B0 shimming workflow.

The pipeline is split between **MATLAB** and **Julia**:

* **MATLAB** loads the raw scanner data and performs any vendor-specific reconstruction steps.
* **Julia** estimates an unwrapped and regularized B0 field map using modern open-source packages, including **ROMEO** for phase unwrapping and **MRIFieldmaps.jl** for regularized field map estimation. The results can be saved as both MATLAB (`.mat`) and NIfTI (`.nii`/`.nii.gz`) files.

## Workflow

```text
Raw scanner data
        │
        ▼
MATLAB
  - Load raw data
  - Reconstruct k-space/images
  - Save reconstruction_input.mat
        │
        ▼
Julia
  - Load reconstruction_input.mat
  - Estimate unwrapped B0 field map
  - Regularize field map
  - Save fieldmap.mat
  - (Optional) Write NIfTI images
```

## MATLAB

The MATLAB code produces a `reconstruction_input.mat` file containing the reconstructed data and acquisition metadata required for field map estimation.

## Julia

The `julia/` directory contains the field map estimation code and its own Julia environment (`Project.toml` and `Manifest.toml`).

To run the Julia portion:

```bash
cd julia
julia
```

Within the Julia REPL:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("estimateFieldMap.jl")
```

For a description of the required input file, output files, package setup, and usage, see:

```text
julia/README.md
```


