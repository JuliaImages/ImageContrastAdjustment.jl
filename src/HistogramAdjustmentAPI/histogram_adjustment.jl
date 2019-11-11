# usage example for package developer:
#
#     import HistogramAdjustmentAPI: AbstractHistogramAdjustmentAlgorithm,
#                             adjust_histogram, adjust_histogram!

"""
    AbstractHistogramAdjustmentAlgorithm <: AbstractImageFilter

The root type for `ImageContrastAdjustment` package.

Any concrete histogram adjustment algorithm shall subtype it to support
[`adjust_histogram`](@ref) and [`adjust_histogram!`](@ref) APIs.

# Examples

All histogram adjustment algorithms in ImageContrastAdjustment are called in the
following pattern:

```julia
# first generate an algorithm instance
f = LinearStretching()

# then pass the algorithm to `adjust_histogram`
img_adjusted = adjust_histogram(img, f)

# or use in-place version `adjust_histogram!`
img_adjusted = similar(img)
adjust_histogram!(img_adjusted, img, f)
```

Some algorithms also receive additional information as an argument,
e.g., `nbins` of `Equalization`.

```julia
# you can explicit specify the parameters
f = Equalization(nbins = 32)
```

For more examples, please check [`adjust_histogram`](@ref),
[`adjust_histogram!`](@ref) and concrete algorithms.
"""
abstract type AbstractHistogramAdjustmentAlgorithm <: AbstractImageFilter end

adjust_histogram!(out::Union{GenericGrayImage, AbstractArray{<:Color3}},
          img,
          f::AbstractHistogramAdjustmentAlgorithm,
          args...; kwargs...) =
    f(out, img, args...; kwargs...)

# TODO: Relax this to all color types
function adjust_histogram!(img::Union{GenericGrayImage, AbstractArray{<:Color3}},
                   f::AbstractHistogramAdjustmentAlgorithm,
                   args...; kwargs...)
    tmp = copy(img)
    f(img, tmp, args...; kwargs...)
    return img
end

function adjust_histogram(::Type{T},
                  img,
                  f::AbstractHistogramAdjustmentAlgorithm,
                  args...; kwargs...) where T
    out = similar(Array{T}, axes(img))
    adjust_histogram!(out, img, f, args...; kwargs...)
    return out
end

adjust_histogram(img::AbstractArray{T},
                 f::AbstractHistogramAdjustmentAlgorithm,
                 args...; kwargs...) where T <: Colorant =
         adjust_histogram(T, img, f, args...; kwargs...)

# Do not promote Number to Gray{<:Number}
adjust_histogram(img::AbstractArray{T},
                 f::AbstractHistogramAdjustmentAlgorithm,
                 args...; kwargs...) where T <: Number =
        adjust_histogram(T, img, f, args...; kwargs...)


# Handle instance where the input is a sequence of images.
adjust_histogram!(out_sequence::Vector{T},
          img_sequence,
          f::AbstractHistogramAdjustmentAlgorithm,
          args...; kwargs...) where T <: Union{GenericGrayImage, AbstractArray{<:Color3}} =
    f(out_sequence, img_sequence, args...; kwargs...)

# TODO: Relax this to all color types
function adjust_histogram!(img_sequence::Vector{T},
                   f::AbstractHistogramAdjustmentAlgorithm,
                   args...; kwargs...) where T <: Union{GenericGrayImage, AbstractArray{<:Color3}}
    tmp = copy(img_sequence)
    f(img_sequence, tmp, args...; kwargs...)
    return img_sequence
end

function adjust_histogram(::Type{T},
                  img_sequence::Vector{<:AbstractArray},
                  f::AbstractHistogramAdjustmentAlgorithm,
                  args...; kwargs...) where T
    N  = length(img_sequence)
    out_sequence = [similar(Array{T}, axes(img_sequence[n])) for n = 1:N]
    adjust_histogram!(out_sequence, img_sequence, f, args...; kwargs...)
    return out_sequence
end

adjust_histogram(img_sequence::Vector{<:AbstractArray{T}},
                 f::AbstractHistogramAdjustmentAlgorithm,
                 args...; kwargs...) where T <: Colorant =
         adjust_histogram(T, img_sequence, f, args...; kwargs...)

# Do not promote Number to Gray{<:Number}
adjust_histogram(img_sequence::Vector{<:AbstractArray{T}},
                 f::AbstractHistogramAdjustmentAlgorithm,
                 args...; kwargs...) where T <: Number =
        adjust_histogram(T, img_sequence, f, args...; kwargs...)

### Docstrings

"""
    adjust_histogram!([out,] img, f::AbstractHistogramAdjustmentAlgorithm, args...; kwargs...)

Adjust histogram of `img` using algorithm `f`.

# Output

If `out` is specified, it will be changed in place. Otherwise `img` will be changed in place.

# Examples

Just simply pass an algorithm to `adjust_histogram!`:

```julia
img_adjusted = similar(img)
adjust_histogram!(img_adjusted, img, f)
```

For cases you just want to change `img` in place, you don't necessarily need to manually
allocate `img_adjusted`; just use the convenient method:

```julia
adjust_histogram!(img, f)
```

See also: [`adjust_histogram`](@ref)
"""
adjust_histogram!

"""
    adjust_histogram([T::Type,] img, f::AbstractHistogramAdjustmentAlgorithm, args...; kwargs...)

Adjust histogram of `img` using algorithm `f`.

# Output

The return image `img_adjusted` is an `Array{T}`.

If `T` is not specified, then it's inferred.
# Examples

Just simply pass the input image and algorithm to `adjust_histogram`

```julia
img_adjusted = adjust_histogram(img, f)
```

This reads as "`adjust_histogram` of image `img` using algorithm `f`".

You can also explicitly specify the return type:

```julia
img_adjusted_float32 = adjust_histogram(Gray{Float32}, img, f)
```

See also [`adjust_histogram!`](@ref) for in-place histogram adjustment.
"""
adjust_histogram
