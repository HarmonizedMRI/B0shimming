# B0 field map reconstruction

This folder contains the reconstruction pipeline for the vendor-neutral B0 shimming workflow.

The pipeline is split between **MATLAB** and **Julia**:

- **MATLAB** loads the raw scanner data, sorts the acquired k-space data, and prepares the input for field map estimation. Optionally, it can also calculate ESPIRiT receive-coil sensitivity maps using BART.
- **Julia** reconstructs the images and estimates an unwrapped and regularized B0 field map using **ROMEO** for phase unwrapping and **MRIFieldmaps.jl** for regularized field map estimation.

## Workflow

```text
Raw scanner data
        │
        ▼
MATLAB
  - Load GE ScanArchive data
  - Sort acquired k-space
  - Save reconstruction_input.mat
  - Optionally calculate ESPIRiT sensitivity maps (sens.mat)
        │
        ▼
Julia
  - Reconstruct images
  - Estimate unwrapped B0 field map
  - Regularize field map
  - Save fieldmap_results.mat
  - Write fieldmap and magnitude NIfTI images
```

## Run the complete reconstruction

The complete MATLAB + Julia pipeline can be run from the `reconstruction/` directory:

```bash
./run_reconstruct_fieldmap.sh data.h5 output
```

This creates:

```text
output/
├── reconstruction_input.mat
├── fieldmap_results.mat
├── fieldmap_hz.nii.gz
└── magnitude_echo1.nii.gz
```

Sensitivity maps can optionally be calculated using BART ESPIRiT:

```bash
./run_reconstruct_fieldmap.sh --sens data.h5 output
```

which additionally creates:

```text
output/
└── sens.mat
```

BART must be installed and its MATLAB interface available on the MATLAB path when using this option.

## MATLAB

The MATLAB code is in `matlab/`.

`reconstructB0.m` accepts either a GE ScanArchive filename or a cell array of raw acquisitions such as that returned by `pge2.utils.loaddata()`.

For example:

```matlab
reconstructB0('data.h5', 'reconstruction_input.mat');
```

Sensitivity maps can optionally be calculated using BART ESPIRiT:

```matlab
reconstructB0( ...
    'data.h5', ...
    'reconstruction_input.mat', ...
    'CalculateSens', true);
```

The resulting `reconstruction_input.mat` contains the k-space data, echo times, and acquisition geometry needed by the Julia reconstruction.

## Julia

The `julia/` directory contains the field map estimation code and its own Julia environment (`Project.toml` and `Manifest.toml`).

The Julia processing can also be run separately for development or debugging. See:

```text
julia/README.md
```

for environment setup, interactive Julia usage, input/output details, and NIfTI writing.