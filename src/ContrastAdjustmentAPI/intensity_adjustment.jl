# usage example for package developer:
#
#     import ContrastAdjustmentAPI: AbstractIntensityAdjustmentAlgorithm,
#                             adjust_intensity, adjust_intensity!

"""
    AbstractIntensityAdjustmentAlgorithm <: AbstractImageFilter

A root type for `ImageContrastAdjustment` package that relates to algorithms
that manipulate contrast without operating on intensity histograms.

Any concrete intensity adjustment algorithm shall subtype it to support
[`adjust_intensity`](@ref) and [`adjust_intensity!`](@ref) APIs.

# Examples

All intensity adjustment algorithms in ImageContrastAdjustment are called in the
following pattern:

```julia
# first generate an algorithm instance
f = LinearStretching()

# then pass the algorithm to `adjust_intensity`
img_adjusted = adjust_intensity(img, f)

# or use in-place version `adjust_intensity!`
img_adjusted = similar(img)
adjust_intensity!(img_adjusted, img, f)
```

Some algorithms also receive additional information as an argument,
e.g., `gamma` of `GammaCorrection`.

```julia
# you can explicit specify the parameters
f = GammaCorrection(gamma = 1.4)
```

For more examples, please check [`adjust_intensity`](@ref),
[`adjust_intensity!`](@ref) and concrete algorithms.
"""
abstract type AbstractIntensityAdjustmentAlgorithm <: AbstractImageFilter end

adjust_intensity!(out::Union{GenericGrayImage, AbstractArray{<:Color3}},
          img,
          f::AbstractIntensityAdjustmentAlgorithm,
          args...; kwargs...) =
    f(out, img, args...; kwargs...)

# TODO: Relax this to all color types
function adjust_intensity!(img::Union{GenericGrayImage, AbstractArray{<:Color3}},
                   f::AbstractIntensityAdjustmentAlgorithm,
                   args...; kwargs...)
    tmp = copy(img)
    f(img, tmp, args...; kwargs...)
    return img
end

function adjust_intensity(::Type{T},
                  img,
                  f::AbstractIntensityAdjustmentAlgorithm,
                  args...; kwargs...) where T
    out = similar(Array{T}, axes(img))
    adjust_intensity!(out, img, f, args...; kwargs...)
    return out
end

adjust_intensity(img::AbstractArray{T},
                 f::AbstractIntensityAdjustmentAlgorithm,
                 args...; kwargs...) where T <: Colorant =
         adjust_intensity(T, img, f, args...; kwargs...)

# Do not promote Number to Gray{<:Number}
adjust_intensity(img::AbstractArray{T},
                 f::AbstractIntensityAdjustmentAlgorithm,
                 args...; kwargs...) where T <: Number =
        adjust_intensity(T, img, f, args...; kwargs...)


# Handle instance where the input is a sequence of images.
adjust_intensity!(out_sequence::Vector{T},
          img_sequence,
          f::AbstractIntensityAdjustmentAlgorithm,
          args...; kwargs...) where T <: Union{GenericGrayImage, AbstractArray{<:Color3}} =
    f(out_sequence, img_sequence, args...; kwargs...)

# TODO: Relax this to all color types
function adjust_intensity!(img_sequence::Vector{T},
                   f::AbstractIntensityAdjustmentAlgorithm,
                   args...; kwargs...) where T <: Union{GenericGrayImage, AbstractArray{<:Color3}}
    tmp = copy(img_sequence)
    f(img_sequence, tmp, args...; kwargs...)
    return img_sequence
end

function adjust_intensity(::Type{T},
                  img_sequence::Vector{<:AbstractArray},
                  f::AbstractIntensityAdjustmentAlgorithm,
                  args...; kwargs...) where T
    N  = length(img_sequence)
    out_sequence = [similar(Array{T}, axes(img_sequence[n])) for n = 1:N]
    adjust_intensity!(out_sequence, img_sequence, f, args...; kwargs...)
    return out_sequence
end

adjust_intensity(img_sequence::Vector{<:AbstractArray{T}},
                 f::AbstractIntensityAdjustmentAlgorithm,
                 args...; kwargs...) where T <: Colorant =
         adjust_intensity(T, img_sequence, f, args...; kwargs...)

# Do not promote Number to Gray{<:Number}
adjust_intensity(img_sequence::Vector{<:AbstractArray{T}},
                 f::AbstractIntensityAdjustmentAlgorithm,
                 args...; kwargs...) where T <: Number =
        adjust_intensity(T, img_sequence, f, args...; kwargs...)

### Docstrings

"""
    adjust_intensity!([out,] img, f::AbstractIntensityAdjustmentAlgorithm, args...; kwargs...)

Adjust intensity of `img` using algorithm `f`.

# Output

If `out` is specified, it will be changed in place. Otherwise `img` will be changed in place.

# Examples

Just simply pass an algorithm to `adjust_intensity!`:

```julia
img_adjusted = similar(img)
adjust_intensity!(img_adjusted, img, f)
```

For cases you just want to change `img` in place, you don't necessarily need to manually
allocate `img_adjusted`; just use the convenient method:

```julia
adjust_intensity!(img, f)
```

See also: [`adjust_intensity`](@ref)
"""
adjust_intensity!

"""
    adjust_intensity([T::Type,] img, f::AbstractIntensityAdjustmentAlgorithm, args...; kwargs...)

Adjust intensity of `img` using algorithm `f`.

# Output

The return image `img_adjusted` is an `Array{T}`.

If `T` is not specified, then it's inferred.
# Examples

Just simply pass the input image and algorithm to `adjust_intensity`

```julia
img_adjusted = adjust_intensity(img, f)
```

This reads as "`adjust_intensity` of image `img` using algorithm `f`".

You can also explicitly specify the return type:

```julia
img_adjusted_float32 = adjust_intensity(Gray{Float32}, img, f)
```

See also [`adjust_intensity!`](@ref) for in-place intensity adjustment.
"""
adjust_intensity
