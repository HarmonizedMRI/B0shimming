#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

include(joinpath(@__DIR__, "src", "B0Reconstruction.jl"))
using .B0Reconstruction

function usage()
    println("""
    Usage:
        julia --project=. estimateFieldMap.jl INPUT.mat OUTPUT.mat

    INPUT.mat must contain:
        kspace      Complex array [nx, nEcho*ny, nz, ncoil]
        echo_times  Echo times in seconds [nEcho]

    Optional fields:
        mask        Boolean or numeric array [nx, ny, nz]
        mask_threshold
                    Relative magnitude threshold used when mask is absent
                    (default: 0.40)

    OUTPUT.mat contains:
        fieldmap_hz
        unwrapped_fieldmap_hz
        initial_fieldmap_hz
        magnitude
        mask
        echo_times
    """)
end

if length(ARGS) != 2
    usage()
    exit(1)
end

input_filename, output_filename = ARGS

result = estimate_fieldmap_file(input_filename, output_filename)

println("Wrote regularized B0 field map to:")
println("  ", abspath(output_filename))
println("Field-map size: ", size(result.fieldmap_hz))
