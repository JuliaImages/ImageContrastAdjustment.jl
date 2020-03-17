
"""
```
    LinearStretching <: AbstractHistogramAdjustmentAlgorithm
    LinearStretching(; minval = 0, maxval = 1, [src_minval], [src_maxval])

    adjust_histogram([T,] img, f::LinearStretching)
    adjust_histogram!([out,] img, f::LinearStretching)
```

Returns an image where the range of the intensities spans the interval [`minval`, `maxval`].

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

## Choices for `minval` and `maxval`

If minval and maxval are specified then intensities are mapped to the range
[`minval`, `maxval`]. The default values are 0 and 1.

## Choices for `src_minval` and `src_maxval`

`src_minval` and `src_maxval` specifies the intensity range of input image. By default,
the values are `extrema(img)` (finite). If custom values are provided, the output
intensity value will be clamped to range `(minval, maxval)` if it exceeds that.

# Example

```julia
using ImageContrastAdjustment, TestImages

img = testimage("mandril_gray")
imgo = adjust_histogram(img, LinearStretching(minval = 0, maxval = 1))

```

# References
1. W. Burger and M. J. Burge. *Digital Image Processing*. Texts in Computer Science, 2016. [doi:10.1007/978-1-4471-6684-9](https://doi.org/10.1007/978-1-4471-6684-9)

"""
@with_kw struct LinearStretching{T₁ <: Union{Real,AbstractGray},
                                 T₂ <: Union{Real,AbstractGray},
                                 T₃ <: Union{Nothing, Real,AbstractGray},
                                 T₄ <: Union{Nothing, Real,AbstractGray}} <: AbstractHistogramAdjustmentAlgorithm
    minval::T₁     = 0.0
    maxval::T₂     = 1.0
    src_minval::T₃ = nothing
    src_maxval::T₄ = nothing
end

function (f::LinearStretching)(out::GenericGrayImage, img::GenericGrayImage)
    img_min, img_max = minfinite(img), maxfinite(img)
    src_minval = isnothing(f.src_minval) ? img_min : f.src_minval
    src_maxval = isnothing(f.src_maxval) ? img_max : f.src_maxval
    T = eltype(out)

    # r * x - o is equivalent to (x-A) * ((b-a)/(B-A)) + a
    # precalculate these and make inner loop contains only multiplication and addition
    # to get better performance
    r = (f.maxval - f.minval) / (src_maxval - src_minval)
    o = (src_minval*f.maxval - src_maxval*f.minval) / (src_maxval - src_minval)

    if 1 ≈ r && 0 ≈ o
        # when image intensity is already adjusted, there's no need to do it again
        # it's a trivial but common case in practice
        out .= img
        return out
    end

    # In most cases, we don't need to clamp the output
    # this is only used when user specifies custom `(src_minval, src_maxval)`
    out_minval = r * img_min - o
    out_maxval = r * img_max - o
    do_clamp = (out_minval < f.minval) || (out_maxval > f.maxval)

    @inbounds @simd for p in eachindex(img)
        val = img[p]
        if isnan(val)
            out[p] = val
        else
            newval = r * val - o
            do_clamp && (newval = clamp(newval, f.minval, f.maxval))
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
