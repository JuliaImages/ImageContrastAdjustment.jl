
"""
```
    LinearStretching <: AbstractHistogramAdjustmentAlgorithm
    LinearStretching(; [src_minval], [src_maxval],
                       dst_minval=0, dst_maxval=1,
                       no_clamp=false)

    LinearStretching((src_minval, src_maxval) => (dst_minval, dst_maxval))
    LinearStretching((src_minval, src_maxval) => nothing)
    LinearStretching(nothing => (dst_minval, dst_maxval))

    adjust_histogram([T,] img, f::LinearStretching)
    adjust_histogram!([out,] img, f::LinearStretching)
```

Returns an image where the range of the intensities spans the interval [`dst_minval`, `dst_maxval`].

# Details


Linear stretching (also called *normalization*) is a contrast enhancing
transformation that is used to modify the dynamic range of the image. In
particular, suppose that the input image has gray values in the range [A,B] and
one wishes to change the dynamic range to [a,b] using a linear mapping, then the
necessary transformation is given by the relation

```math
f(x) = (x-A) \\frac{b-a}{B-A} + a.
```

# Options

Various options for the parameters of the `adjust_histogram` and
`LinearStretching` type  are described in more detail below.

## Choices for `img`

The function can handle a variety of input types. The returned
image depends on the input type.

For colored images, the input is converted to the
[YIQ](https://en.wikipedia.org/wiki/YIQ)  type and the intensities of the Y
channel are stretched to the specified range. The modified Y channel is then
combined with the I and Q channels and the resulting image converted to the same
type as the input.

## Choices for `dst_minval` and `dst_maxval`

If destination value range `dst_minval` and `dst_maxval` are specified then intensities are
mapped to the range [`dst_minval`, `dst_maxval`]. The default values are 0 and 1.

## Choices for `src_minval` and `src_maxval`

The source value range `src_minval` and `src_maxval` specifies the intensity range of input
image. By default, the values are `extrema(img)` (finite). If custom values are provided,
the output intensity value will be clamped to range `[dst_minval, dst_maxval]` if it exceeds that.

## `no_clamp`

Setting `no_clamp=true` to disable the automatic clamp even if the output intensity value
exceeds the range `[dst_minval, dst_maxval]`. Note that a clamp is still applied for types
that has limited value range, for example, if the input eltype is `N0f8`, then the output will
be clamped to `[0.0N0f8, 1.0N0f8]` even if `no_clamp==true`.

# Example

```julia
using ImageContrastAdjustment, TestImages

img = testimage("mandril_gray")
# Stretches the contrast in `img` so that it spans the unit interval.
imgo = adjust_histogram(img, LinearStretching(dst_minval = 0, dst_maxval = 1))
```

For convenience, Constructing a `LinearStretching` object using `Pair` is also supported

```julia
# these two constructors are equivalent
LinearStretching(src_minval=0.1, src_maxval=0.9, dst_minval=0.05, dst_maxval=0.95)
LinearStretching((0.1, 0.9) => (0.05, 0.95))

# replace the part with `nothing` to use default values, e.g.,
# specify only destination value range
LinearStretching(nothing => (0.05, 0.95))
# specify only source value range and use default destination value range, i.e., (0, 1)
LinearStretching((0.1, 0.9) => nothing)
```

# References
1. W. Burger and M. J. Burge. *Digital Image Processing*. Texts in Computer Science, 2016. [doi:10.1007/978-1-4471-6684-9](https://doi.org/10.1007/978-1-4471-6684-9)

"""
@with_kw struct LinearStretching{T} <: AbstractHistogramAdjustmentAlgorithm
    src_minval::T = nothing
    src_maxval::T = nothing
    dst_minval::T = 0.0f0
    dst_maxval::T = 1.0f0
    minval::T = nothing
    maxval::T = nothing
    no_clamp::Bool = false
    function LinearStretching(src_minval::T1,
                              src_maxval::T2,
                              dst_minval::T3,
                              dst_maxval::T4,
                              minval::T5=nothing,
                              maxval::T6=nothing,
                              no_clamp::Bool=false) where {T1 <: Union{Nothing,Real,AbstractGray},
                                                 T2 <: Union{Nothing,Real,AbstractGray},
                                                 T3 <: Union{Nothing,Real,AbstractGray},
                                                 T4 <: Union{Nothing,Real,AbstractGray},
                                                 T5 <: Union{Nothing,Real,AbstractGray},
                                                 T6 <: Union{Nothing,Real,AbstractGray}}
        # in order to deprecate old fields we have to introduce new fields if we still want to use @with_kw
        # https://github.com/JuliaImages/ImageContrastAdjustment.jl/pull/28#discussion_r395751301
        if !isnothing(minval)
            dst_minval = minval
            Base.depwarn("deprecated: use `dst_minval` for keyword `minval`", :LinearStretching)
        end
        if !isnothing(maxval)
            dst_maxval = maxval
            Base.depwarn("deprecated: use `dst_maxval` for keyword `maxval`", :LinearStretching)
        end

        dst_minval <= dst_maxval || throw(ArgumentError("dst_minval $dst_minval should be less than dst_maxval $dst_maxval"))
        if !(isnothing(src_minval) || isnothing(src_maxval))
            src_minval <= src_maxval || throw(ArgumentError("src_minval $src_minval should be less than src_maxval $src_maxval"))
        end
        T = promote_type(T1, T2, T3, T4, T5, T6)
        new{T}(convert(T, src_minval), convert(T, src_maxval),
               convert(T, dst_minval), convert(T, dst_maxval),
               convert(T, dst_minval), convert(T, dst_maxval),
               no_clamp)
    end
end
function LinearStretching(rangemap::Pair{Tuple{T1, T2}, Tuple{T3, T4}}; no_clamp=false) where {T1, T2, T3, T4}
    LinearStretching(rangemap.first..., rangemap.second..., nothing, nothing, no_clamp)
end
function LinearStretching(rangemap::Pair{Nothing, Tuple{T3, T4}}; no_clamp=false) where {T3, T4}
    LinearStretching(nothing, nothing, rangemap.second..., nothing, nothing, no_clamp)
end
function LinearStretching(rangemap::Pair{Tuple{T1, T2}, Nothing}; no_clamp=false) where {T1, T2}
    LinearStretching(src_minval=rangemap.first[1], src_maxval=rangemap.first[2], no_clamp=no_clamp)
end

function (f::LinearStretching)(out::GenericGrayImage, img::GenericGrayImage)
    T = eltype(out)
    FT = eltype(floattype(T))
    img_min, img_max = minimum_finite(img), maximum_finite(img)
    # explicit annotation is needed because the ?: line mixes three value types:
    # Nothing, T, and typeof(f.src_minval)
    src_minval::FT = isnothing(f.src_minval) ? img_min : f.src_minval
    src_maxval::FT = isnothing(f.src_maxval) ? img_max : f.src_maxval
    dst_minval::FT = f.dst_minval
    dst_maxval::FT = f.dst_maxval

    # the kernel operation `r * x - o` is equivalent to `(x-A) * ((b-a)/(B-A)) + a`
    # precalculate these and make inner loop contains only multiplication and addition
    # to get better performance
    r = (dst_maxval - dst_minval) / (src_maxval - src_minval)
    o = (src_minval*dst_maxval - src_maxval*dst_minval) / (src_maxval - src_minval)

    if 1 ≈ r && 0 ≈ o
        # when image intensity is already adjusted, there's no need to do it again
        # it's a trivial but common case in practice
        out === img || (out .= img)
        return out
    end

    # In most cases, we don't need to clamp the output
    # this is only used when user specifies custom parameters
    out_minval = r * img_min - o
    out_maxval = r * img_max - o
    do_clamp = f.no_clamp ? false : (out_minval < dst_minval) || (out_maxval > dst_maxval)
    # re-adjust clamp option to avoid ArgumentError
    if !do_clamp && (out_minval < typemin(T) || out_maxval > typemax(T))
        do_clamp = true
        dst_minval = convert(typeof(dst_minval), typemin(eltype(T)))
        dst_maxval = convert(typeof(dst_maxval), typemax(eltype(T)))
    else
        do_clamp = do_clamp || out_minval < typemin(T) || out_maxval > typemax(T)
        dst_minval = max(dst_minval, convert(typeof(dst_minval), typemin(eltype(T))))
        dst_maxval = min(dst_maxval, convert(typeof(dst_maxval), typemax(eltype(T))))
    end

    # although slightly faster, clamp before linear stretching has some rounding issue, e.g,
    # it's likely to get results like -9.313226f-10, which can't be directly converted to N0f8
    # thus the line below is commented out and not used
    # do_clamp && (img = clamp.(img, src_minval, src_maxval))

    # tweak the performance of FixedPoint by fusing operations into one broadcast
    # for Float32 the fallback implementation is faster
    if eltype(T) <: FixedPoint
        # ?: is faster than if-else
        @. out = do_clamp ? clamp(r * img - o, dst_minval, dst_maxval) : r * img - o
        return out
    end

    # fallback implementation
    @inbounds @simd for p in eachindex(img)
        val = img[p]
        if isnan(val)
            out[p] = val
        else
            newval = r * val - o
            do_clamp && (newval = clamp(newval, dst_minval, dst_maxval))
            out[p] = T <: Integer ? round(Int, newval) : newval
        end
    end
    out
end

function (f::LinearStretching)(out::AbstractArray{<:Color3}, img::AbstractArray{<:Color3})
    T = eltype(out)
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    adjust_histogram!(view(yiq_view,1,:,:), f)
    out .= convert.(T, yiq)
end

(f::LinearStretching)(out::GenericGrayImage, img::AbstractArray{<:Color3}) =
    f(out, of_eltype(Gray, img))
