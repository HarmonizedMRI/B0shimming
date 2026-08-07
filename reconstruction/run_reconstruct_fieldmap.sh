#!/usr/bin/env bash

set -euo pipefail

# -------------------------------------------------------------------------
# Usage
# -------------------------------------------------------------------------

if [ "$#" -lt 2 ]; then
    echo "Usage:"
    echo "  $0 <ScanArchive.h5> <output_directory>"
    echo
    echo "Example:"
    echo "  $0 data.h5 output"
    exit 1
fi

INPUT=$(realpath "$1")
OUTPUT_DIR=$(realpath -m "$2")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MATLAB_DIR="${SCRIPT_DIR}/matlab"
JULIA_DIR="${SCRIPT_DIR}/julia"

mkdir -p "${OUTPUT_DIR}"

RECON_INPUT="${OUTPUT_DIR}/reconstruction_input.mat"
FIELDMAP_MAT="${OUTPUT_DIR}/fieldmap_results.mat"
FIELDMAP_NII="${OUTPUT_DIR}/fieldmap_hz.nii.gz"
MAGNITUDE_NII="${OUTPUT_DIR}/magnitude_echo1.nii.gz"


# -------------------------------------------------------------------------
# MATLAB: load ScanArchive and prepare reconstruction input
# -------------------------------------------------------------------------

echo
echo "========================================"
echo "MATLAB reconstruction"
echo "========================================"
echo

matlab -batch "\
    cd('${MATLAB_DIR}'); \
    setup; \
    reconstructB0('${INPUT}', '${RECON_INPUT}');"


# -------------------------------------------------------------------------
# Julia: unwrap and regularize B0 field map
# -------------------------------------------------------------------------

echo
echo "========================================"
echo "Julia field-map estimation"
echo "========================================"
echo

julia --project="${JULIA_DIR}" -e "\
    include(joinpath(\"${JULIA_DIR}\", \"src\", \"B0Reconstruction.jl\")); \
    using .B0Reconstruction; \
    estimate_fieldmap_file( \
        \"${RECON_INPUT}\", \
        \"${FIELDMAP_MAT}\"; \
        fieldmap_nifti=\"${FIELDMAP_NII}\", \
        magnitude_nifti=\"${MAGNITUDE_NII}\" \
    )"


# -------------------------------------------------------------------------
# Done
# -------------------------------------------------------------------------

echo
echo "========================================"
echo "Reconstruction complete"
echo "========================================"
echo
echo "Outputs:"
echo "  ${RECON_INPUT}"
echo "  ${FIELDMAP_MAT}"
echo "  ${FIELDMAP_NII}"
echo "  ${MAGNITUDE_NII}"

