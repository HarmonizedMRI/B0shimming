module B0Reconstruction

using FFTW: fftshift, ifftshift, ifft
using MAT: matread, matwrite
using MRIFieldmaps: b0init, b0map, b0scale
using ROMEO: unwrap

export reconstruct_echoes,
       create_mask,
       estimate_fieldmap,
       estimate_fieldmap_file

"""
    reconstruct_echoes(kspace, n_echoes)

Reconstruct interleaved multi-coil Cartesian k-space.

# Arguments
- `kspace`: complex array with dimensions `[nx, nEcho*ny, nz, ncoil]`.
- `n_echoes`: number of interleaved echoes.

# Returns
A complex array with dimensions `[nx, ny, nz, ncoil, nEcho]`.

The second k-space dimension is assumed to be ordered as

    echo 1, echo 2, ..., echo nEcho,
    echo 1, echo 2, ..., echo nEcho, ...

for successive phase-encoding locations.
"""
function reconstruct_echoes(kspace::AbstractArray, n_echoes::Integer)
    ndims(kspace) == 4 ||
        throw(ArgumentError(
            "kspace must have dimensions [nx, nEcho*ny, nz, ncoil]"
        ))
    n_echoes >= 2 ||
        throw(ArgumentError("At least two echoes are required"))

    nx, ny_interleaved, nz, ncoil = size(kspace)
    ny_interleaved % n_echoes == 0 ||
        throw(ArgumentError(
            "The second k-space dimension ($ny_interleaved) is not " *
            "divisible by n_echoes ($n_echoes)"
        ))

    ny = ny_interleaved ÷ n_echoes
    images = Array{ComplexF32}(undef, nx, ny, nz, ncoil, n_echoes)

    # Perform a 3D inverse FFT independently for every coil and echo.
    for echo_index in 1:n_echoes
        echo_kspace = ComplexF32.(
            @view kspace[:, echo_index:n_echoes:end, :, :]
        )

        for coil_index in 1:ncoil
            coil_kspace = @view echo_kspace[:, :, :, coil_index]
            images[:, :, :, coil_index, echo_index] =
                fftshift(
                    ifft(ifftshift(coil_kspace, (1, 2, 3)), (1, 2, 3)),
                    (1, 2, 3),
                )
        end
    end

    return images
end

"""
    create_mask(magnitude; threshold=0.40)

Create a simple mask by thresholding a magnitude image relative to its maximum.
"""
function create_mask(
    magnitude::AbstractArray;
    threshold::Real = 0.40,
)
    0 <= threshold <= 1 ||
        throw(ArgumentError("threshold must be between 0 and 1"))

    max_magnitude = maximum(magnitude)
    max_magnitude > 0 ||
        throw(ArgumentError("magnitude image is identically zero"))

    return BitArray(magnitude .>= threshold * max_magnitude)
end

"""
    estimate_fieldmap(images, echo_times; mask=nothing,
                      mask_threshold=0.40, l2b=-3, precon=:diag)

Estimate an unwrapped and regularized B0 field map.

# Arguments
- `images`: complex multi-coil images with dimensions
  `[nx, ny, nz, ncoil, nEcho]`.
- `echo_times`: echo times in seconds.
- `mask`: optional 3D processing mask.
- `mask_threshold`: relative magnitude threshold when `mask` is omitted.
- `l2b`: regularization parameter passed to `MRIFieldmaps.b0map`.
- `precon`: preconditioner passed to `MRIFieldmaps.b0map`.

# Returns
A named tuple containing the regularized field map, unwrapped initial field
map, wrapped initial field map, magnitude image, mask, and echo times.
All field maps are in Hz.
"""
function estimate_fieldmap(
    images::AbstractArray,
    echo_times::AbstractVector;
    mask = nothing,
    mask_threshold::Real = 0.40,
    l2b::Real = -3,
    precon::Symbol = :diag,
)
    ndims(images) == 5 ||
        throw(ArgumentError(
            "images must have dimensions [nx, ny, nz, ncoil, nEcho]"
        ))

    n_echoes = size(images, 5)
    length(echo_times) == n_echoes ||
        throw(DimensionMismatch(
            "length(echo_times) must equal size(images, 5)"
        ))
    n_echoes == 2 ||
        throw(ArgumentError(
            "This initial implementation currently requires exactly two echoes"
        ))

    te = Float32.(vec(echo_times))
    all(isfinite, te) ||
        throw(ArgumentError("echo_times must be finite"))
    issorted(te) ||
        throw(ArgumentError("echo_times must be increasing"))

    delta_te = te[2] - te[1]
    delta_te > 0 ||
        throw(ArgumentError("Echo-time difference must be positive"))

    # Preserve the coil combination used in the existing shimtool.jl:
    # complex summation across receive coils.
    combined = sum(ComplexF32.(images); dims=4)

    magnitude = dropdims(
        sqrt.(sum(abs2, @view(images[:, :, :, :, 1]); dims=4));
        dims=4,
    )

    processing_mask = if isnothing(mask)
        create_mask(magnitude; threshold=mask_threshold)
    else
        size(mask) == size(magnitude) ||
            throw(DimensionMismatch(
                "mask size $(size(mask)) does not match image size " *
                "$(size(magnitude))"
            ))
        BitArray(mask .!= 0)
    end

    any(processing_mask) ||
        throw(ArgumentError("processing mask is empty"))

    # MRIFieldmaps scaling and wrapped initialization.
    scaled_combined, scale = b0scale(combined, te)
    scaled_images = ComplexF32.(images ./ scale)

    initial_result = b0init(scaled_images, te)

    initial_hz = if ndims(initial_result) == 4 && size(initial_result, 4) == 1
        dropdims(initial_result; dims=4)
    elseif ndims(initial_result) == 3
        initial_result
    else
        throw(DimensionMismatch(
            "Unexpected output size from b0init: $(size(initial_result))"
        ))
    end

    initial_hz = Float32.(initial_hz)
    initial_hz .*= processing_mask

    # Convert the initial map to phase difference, unwrap with ROMEO,
    # then convert back to Hz.
    initial_phase = 2f0 * Float32(pi) * delta_te .* initial_hz
    unwrapped_phase = unwrap(
        initial_phase;
        mag=magnitude,
        mask=processing_mask,
    )
    unwrapped_hz =
        Float32.(unwrapped_phase ./ (2f0 * Float32(pi) * delta_te))
    unwrapped_hz .*= processing_mask

    # Regularized field-map estimation.
    fieldmap_result = b0map(
        scaled_images,
        te;
        smap=nothing,
        l2b=l2b,
        precon=precon,
        finit=unwrapped_hz,
        mask=processing_mask,
    )
    fieldmap_hz = Float32.(fieldmap_result[1])

    return (
        fieldmap_hz=fieldmap_hz,
        unwrapped_fieldmap_hz=unwrapped_hz,
        initial_fieldmap_hz=Float32.(initial_hz),
        magnitude=Float32.(magnitude),
        mask=processing_mask,
        echo_times=te,
    )
end

"""
    estimate_fieldmap_file(input_filename, output_filename; kwargs...)

Load interleaved k-space from a MATLAB file, reconstruct the echoes, estimate
the field map, and write the results to another MATLAB file.
"""
function estimate_fieldmap_file(
    input_filename::AbstractString,
    output_filename::AbstractString;
    l2b::Real = -3,
    precon::Symbol = :diag,
)
    isfile(input_filename) ||
        throw(ArgumentError("Input file does not exist: $input_filename"))

    input = matread(input_filename)

    haskey(input, "kspace") ||
        throw(KeyError(
            "Input MAT file must contain a variable named \"kspace\""
        ))
    haskey(input, "echo_times") ||
        throw(KeyError(
            "Input MAT file must contain a variable named \"echo_times\""
        ))

    kspace = input["kspace"]
    echo_times = Float32.(vec(input["echo_times"]))
    n_echoes = length(echo_times)

    mask_threshold = if haskey(input, "mask_threshold")
        Float64(input["mask_threshold"])
    else
        0.40
    end

    mask = get(input, "mask", nothing)

    images = reconstruct_echoes(kspace, n_echoes)
    result = estimate_fieldmap(
        images,
        echo_times;
        mask=mask,
        mask_threshold=mask_threshold,
        l2b=l2b,
        precon=precon,
    )

    output = Dict(
        "fieldmap_hz" => result.fieldmap_hz,
        "unwrapped_fieldmap_hz" => result.unwrapped_fieldmap_hz,
        "initial_fieldmap_hz" => result.initial_fieldmap_hz,
        "magnitude" => result.magnitude,
        "mask" => UInt8.(result.mask),
        "echo_times" => result.echo_times,
    )
    matwrite(output_filename, output; compress=true)

    return result
end

end
